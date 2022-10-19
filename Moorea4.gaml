

model Moorea

global {
	
	//-------paramètres du modele---------------------------------------------------------------------------------------------------------------------------
	
	//nous dit si l'on souhaite ou non utiliser les donnees clippées ou complete
	bool clipped_data <- false;
	
	// AMP
	bool AMP <- true;
	
	//est-ce que l'on souhaite activer le pre-calcule et la sauvegarde des distances (à ne faire qu'une fois)
	bool pre_calculer_distances <- false;
	string chemin_fichier_distances <- "../includes/distances.csv";
	

	//a partir de combien de cellules d'eau contigues les conserve-t-on dans la liste des cellules ou l'on peut pecher
	int seuil_filtrage_cellules_isolees <- 20;
	
	
	//duree d'un pas de temps
	float step <- 1.0;  
		
	//duree de la simulation
	float temps_simulation <- 7200.5;
	float temps <- 0.0; // en nombre de jours
	float compteur;
	int temps_annees <- 1;
	
	// Scénarios de gestion	
	float quota <- 10000.0; 
	float aide_financiere <- 0.0; // aide à la mobilité
	bool interdiction_peche_nocturne <- false;	
	
	// COTS
	bool perturbation <- false;
	float debut_perturbation <- 1000.0; // début de la perturbation COTS
	
	// Surveillance	
	float surveillance <- 0.2;
	
	float biomasse_pechee_moyenne <- 5.0;
	float coeff_capture_jour;
	float coeff_capture;
	
	//distance qu'un pecheur peut parcourir en fonction de son equipement
	float distance_bateau <- 10 #km;
	float distance_pirogue <- 3 #km;
	float distance_nage <- 1 #km;
	
	float temps_peche_moyen;
		
	// Outputs
	int conflits_jour_pecheur;  
	int conflits_jour_tourisme;	
	int conflits_nuit;
	int quota_atteint;
	float biomasse_pechee_herb_lagon;  
	float biomasse_pechee_herb_pente_externe;
	float biomasse_pechee_carniv_lagon;  
	float biomasse_pechee_carniv_pente_externe;
	float biomasse_pechee_herb_lagon_r;  
	float biomasse_pechee_herb_pente_externe_r;
	float biomasse_pechee_carniv_lagon_r;  
	float biomasse_pechee_carniv_pente_externe_r;
	float H_moyenne; 	
	float CO_moyenne; 	
	float P_moyenne;
	float H_moyenne_lagon; 	
	float CO_moyenne_lagon; 	
	float P_moyenne_lagon;
	float H_moyenne_pente_externe; 	
	float CO_moyenne_pente_externe; 	
	float P_moyenne_pente_externe;
	float C_moyenne_pente_externe;	
	float C_moyenne_lagon;
	float T_moyenne_pente_externe;	
	float T_moyenne_lagon;
	float H_moyenne_AMP; 
	float CO_moyenne_AMP; 	
	float P_moyenne_AMP;
	float H_moyenne_not_AMP; 	
	float CO_moyenne_not_AMP; 	
	float P_moyenne_not_AMP;	
	float temps_recuperation_corail;
	string sim_id <- "" + debut_perturbation +"_" + quota + "_" + aide_financiere + "_" + seed + "_";
	string output_carte_lagon <- "../Outputs/Statu quo/Sans COTS/Lagon/statu_quo_lagon_";
	string output_carte_commune <- "../Outputs/Statu quo/Sans COTS/Communes/statu_quo_commune_";
	string output_csv_cell <- "../Outputs/Statu quo/Sans COTS/Cellules/statu_quo_cell_";
	string output_csv_pecheur <- "../Outputs/Statu quo/Sans COTS/Pecheurs/statu_quo_pecheur_";
	string output_csv_ecologique <- "../Outputs/Statu quo/Sans COTS/statu_quo_ecologique.csv";
	string resultats_batch <- "../Outputs/Statu quo/Sans COTS/"+ sim_id+"resultats_batch.csv";
	
	
	//---Fichiers de données-------------------------------------------------------------------------------------------------------------------------------	
	file Households_shape <- shape_file('../includes/households.shp'); 
	file District_shape <- shape_file('../includes/district_Commune_Moorea-Maiao.shp');
	file commune_shape <- shape_file('../includes/COMMUNE.shp');
	file pente_externe_file <- file("../includes/Pente_externe_polygone.shp");
	file amp_file <- file('../includes/matrix/amp.shp');		
	file coral_file <- file('../includes/matrix/coral_cover_alignes.asc');
	file clipped_file <- file('../includes/cells_clipped_2.shp');	
	file pass_file <- shape_file ('../includes/Passes/my_lagoon_0.shp');
	string household_csv_file <- clipped_data ? "../includes/matrix/Mobilite_clipped.csv":"../includes/matrix/Mobilite.csv";	
		
	//coefficient de spill over (diffusion sur le voisinage) des poissons
	float coeff_spill_over <- 0.9;
	
	// COTS
	float coeff_destruction_corail_COTS <- 0.00023 * 4;
	
	//-------variables internes au modele---------------------------------------------------------------------------------------------------------------------------
	//cellule actives (celles où il y a des poissons/algues/corails)
	list<cell> active_cells; 
	int Nb_active_cells_pente_externe;
	int Nb_active_cells_lagon;
	int Nb_active_cells_total;
	
	// Critères max
	float max_preference;	
	float max_tourism;
	
	//geometrie du monde
	geometry shape <- envelope(coral_file); // Le raster corail définit l'environnement de travail
	
	//nous donne si on est en mode batch
	bool mode_batch <- false;
	
	//nous donne si on est le jour ou la nuit
	bool nuit <- false;
	
	list<pecheur> pecheurs_de_nuit; 
	list<pecheur> pecheurs_de_jour;	
	float nb_pecheurs_avec_bateau;
	float nb_pecheurs_avec_pirogue; 		
		
	// Paramètres fixes du modèle de Lotka Volterra	
	float C_min <- 1.0;   // 1% de recouvrement
	float T_min <- 1.0;
			
	float alpha_C_global <- 0.001506849 / 100; 	
	float alpha_T_global <- 0.050780822 / 10;
	float gamma_T <- alpha_T_global;
	
	float alpha_H <- 0.0047068493 / 10;
	float gamma_H <- alpha_H;	
	float alpha_CO <- 0.0021643836 / 10;
	float gamma_CO <- alpha_CO;
	float alpha_P <- 0.002226484 / 10;
	float gamma_P <- alpha_P;
	
	
	init {
		if (mode_batch) {
			save "cycle,H_moyenne_lagon,CO_moyenne_lagon,P_moyenne_lagon,H_moyenne_pente_externe,CO_moyenne_pente_externe,P_moyenne_pente_externe,
			biomasse_totale_pechee_herb_jour,biomasse_totale_pechee_herb_nuit,biomasse_totale_pechee_carniv_jour,biomasse_totale_pechee_carniv_nuit,
			nb_spill_over, quota_atteint_jour, quota_atteint_nuit,conflits_jour,conflits_nuit,C_moyenne_lagon,T_moyenne_lagon,
			C_moyenne_pente_externe,T_moyenne_pente_externe" to: resultats_batch ;
		}
		if (clipped_data) {
			loop c over: clipped_file {
				ask cell overlapping c {
					clipped <- true;   // Cellules de la pente externe	
				}				
			}
		}				
		
		loop p over: pass_file {
			ask cell overlapping p {				
				active <- false;   // Cellules de la pente externe						
			}			
		}
		
		active_cells  <- cell where (each.grid_value > -999 and each.active = true  and ((not clipped_data) or (each.clipped = true)));		
		list<cell> to_remove;		
		
		loop g over: pente_externe_file {
			ask cell overlapping g {
				pente_externe <- true;   // Cellules de la pente externe				
			}
		}
		
		matrix preference_matrix <- matrix(csv_file("../includes/matrix/fisher_preference_2.asc", " ", float)) ; 
		matrix tourism_matrix <- matrix(csv_file("../includes/matrix/r_raw_tourism_shipping_alignes.asc", " ", float));
		matrix turf_matrix <- matrix(csv_file("../includes/matrix/turf_cover_alignes_2.asc", " ", float));
		matrix herbivores_matrix <- matrix(csv_file("../includes/matrix/herbivores_biomass_alignes_2.asc", " ", float));
		matrix corallivores_matrix <- matrix(csv_file("../includes/matrix/corallivores_biomass_alignes_2.asc", " ", float));
		matrix carnivores_matrix <- matrix(csv_file("../includes/matrix/carnivores_biomass_alignes_2.asc", " ", float));
						
		ask active_cells parallel: true{
			
			C <- grid_value; // surface en m² (pour une case)																		
			M <- 3.0;	
			
			if (pente_externe) {
				C <- 3 * C;
				T <- (3/4) * C;					
			}			
			else {
				T <- 20.0;	
			}
			carrying_capacity_C_T <- T + C;
			
			if C <= C_min {
				C <- C_min;				
			}	
			
			if T <= T_min {
				T <- T_min;				
			}	
													
			H <- float(herbivores_matrix [grid_x, grid_y]);	// biomasse en kg
			
			if H <= 0 {
				H <- 0.0;				
			}			
			carrying_capacity_H <- 2 * H;
			H_min <- H / 10;
											
			CO <- float(corallivores_matrix [grid_x, grid_y]);	
			if CO <= 0 {
				CO <- 0.0;				
			}	
			carrying_capacity_CO <- 2 * CO;
			CO_min <- CO / 2;
						
												
			P <- float(carnivores_matrix [grid_x, grid_y]);	
			if CO <= 0 {
				P <- 0.0;				
			}	
			carrying_capacity_P <- 2 * P;
			P_min <- P / 10;
						
			alpha_C <- alpha_C_global;	
			alpha_T <- alpha_T_global;
			C_init <- C;		
			T_init <- T;
			H_init <- H;
			CO_init <- CO;
			P_init <- P;	
							
			beta_T_C_2 <- (alpha_T / 2) / C; // coefficients de Lotka Volterra à l'équilibre
			delta_CO_C_2 <- gamma_CO / C;	      		     		
     		beta_C_T_2 <- (alpha_C / 2) / T;
			delta_H_T_2 <- gamma_H / T;			
			beta_T_H_2 <- (alpha_T / 2) / H;				
			beta_C_CO_2 <- (alpha_C / 2) / CO;						
			delta_P_H_CO_2 <- gamma_P / (H + CO);	
			beta_H_P_2 <- (alpha_H / 2) / P;
			beta_CO_P_2 <- alpha_CO / P;					
								
			level_tourism <- float(tourism_matrix[grid_x, grid_y]);  // Niveau de tourisme  // modifié
			if (level_tourism < 0) {
				level_tourism <- 0.0;				
			}	
			preference <- float(preference_matrix[grid_x, grid_y]);  // Fichier de préférence des pêcheurs
			if (preference < 0) {
				preference <- 0.0;
			}					
			
			if (C < 0.0) {
				to_remove<<self;
			}			
			
			color <- #blue;
			
		}
				
		active_cells <- active_cells - to_remove;  
		ask active_cells {is_active <- true;}	
		
		create household_area from: Households_shape;
		
		
		//on supprime les cellules isolees
		do filtrage_cellules_eau_isolees;
        	
		//pre-calcul et sauvegarde dans un fichier des distances entre ports possibles et cellules possibles
		if (pre_calculer_distances) {
			do pre_calcul_distances;
		}
	
		//chargement du fichier des distances 
       	do chargement_fichier_distances;
						
		if not mode_batch {write "Initialisation des cellules finie ";}
		Nb_active_cells_pente_externe <- length(active_cells where (each.pente_externe = true));
		Nb_active_cells_lagon <- length(active_cells where (each.pente_externe = false));
		Nb_active_cells_total <- length(active_cells);
		if not mode_batch {write "Nombre de cellules en pente externe :" + Nb_active_cells_pente_externe;}
		if not mode_batch {write "Nombre de cellules dans le lagon :" + Nb_active_cells_lagon;}
				
		max_tourism <- active_cells max_of each.level_tourism;
		max_preference <- active_cells max_of each.preference;			
		
		ask active_cells {
			tourism_criterion_jour <- 1 - ln (level_tourism / max_tourism + 1) / ln(2); // Critères
			pref_criterion <- ln (preference / max_preference + 1) / ln(2) ;			
		}		
		
					   	
		create amp from: amp_file with: [statut::int(read("STATUS"))] {
			
			if AMP {
				ask active_cells overlapping self {
			 		statut <- myself.statut_AMP;
				}			 
			}
		}
		if not mode_batch {write "creation des amp finie";}
		
		create District from: District_shape with: [IDDist07::int(read('IDDIST07')), ID_commune_district :: int(read('ID_COMMUNE')), bord_mer :: string(read('BORD_MER')),sensitivity_district::float(read('HH_S07_SC')), nb_pecheurs_pro_district::int(read('PRO')), nb_pecheurs_annexes_district::int(read('ANNEXES')), nb_individus_district :: int (read('NB_IND'))] {
			write(IDDist07);
		}
		
		create commune from: commune_shape with: [ID_commune :: int(read('ID_SAU'))];
		
		create household from:csv_file(household_csv_file, true) with:[ID_District::int(get("CONST")), ID_LOG::int(get("LOG")), bateaux::int(get("Bateaux")), pirogues::int(get("Pirogues")), pecheurs_pro::int(get("Pro")), pecheurs_annexes::int(get("Annexes"))]	{
					 				  			
					create pecheur number: pecheurs_pro  {   // On considère les pêcheurs "annexes" qui sont dans le district
					
	  					temps_peche <- 8.0; // en heures
	  					proba_peche <- rnd (2/7, 5/7);	
	  					IDDist07_pecheur <- myself.ID_District;  // numéro de district	  					
	  					ID_LOG_pecheur <- myself.ID_LOG;
	  				
	  					if (myself.bateaux >= 1) {
	  						fishing_radius <- distance_bateau;		  						
	  						nb_pecheurs_avec_bateau <- nb_pecheurs_avec_bateau + 1; 
	  						bateau_pecheur <- true;					
	  					}
	  					else {
	  						if (myself.pirogues >= 1) {
	  							fishing_radius <- distance_pirogue;	  							
	  							nb_pecheurs_avec_pirogue <- nb_pecheurs_avec_pirogue + 1;
	  							pirogue_pecheur <- true;   						
	  						}
	  						else {
	  							fishing_radius <- distance_nage;
	  						}	  					
	  					}	
	  									
	  					do init_pecheur;  	
	  									
	  				}		
	  				  		
	  				create pecheur number: pecheurs_annexes {   // On considère les pêcheurs pro qui sont dans le district	  	  						
	  					temps_peche <- 4.0;	 
	  					proba_peche <- rnd (1/7, 3/7);
	  					IDDist07_pecheur <- myself.ID_District;  
	  					ID_LOG_pecheur <- myself.ID_LOG;
	  				
	  					if (myself.bateaux >= 1) {
	  						fishing_radius <- distance_bateau;		  						
	  						nb_pecheurs_avec_bateau <- nb_pecheurs_avec_bateau + 1;  
	  						bateau_pecheur <- true;						
	  					}
	  					else {
	  						if (myself.pirogues >= 1) {
	  							fishing_radius <- distance_pirogue;	  							
	  							nb_pecheurs_avec_pirogue <- nb_pecheurs_avec_pirogue + 1;
	  							pirogue_pecheur <- true;     						
	  						}
	  						else {
	  							fishing_radius <- distance_nage;
	  						}	  					
	  					}	  					
	  						  					
	  					do init_pecheur;	               	    				  				
	  				}
	  				if not mode_batch {write "Creation des pecheurs : " + length(pecheur);}
		}			
				
		pecheurs_de_nuit <- pecheur where (each.peche_nuit = true and not empty(each.fishing_area));
	  	if not mode_batch {write("Nb pecheurs de nuit : " + length(pecheurs_de_nuit));}  			
	  	pecheurs_de_jour <- pecheur where (each.peche_nuit = false and not empty(each.fishing_area));
	  	if not mode_batch {write("Nb pecheurs de jour : " + length(pecheurs_de_jour));}
	  	if not mode_batch {	write("Nb pêcheurs avec bateau : " + nb_pecheurs_avec_bateau);}
	  	if not mode_batch {write("Nb pêcheurs avec pirogue : " + nb_pecheurs_avec_pirogue);	}	  	
	  	
	  	ask active_cells {
	  		if (length(pecheur) != 0) {
	  			beta_H_PE <- (alpha_H / 2) / (length(pecheur) / Nb_active_cells_total);
	  			beta_P_PE <-  alpha_P / (length(pecheur) / Nb_active_cells_total);	  			
	  		}	  		
	  		else {
	  			beta_H_PE <- 0.0;
	  			beta_P_PE <-  0.0;	  
	  		}
	  	}
	  	
	  	H_moyenne <- mean(active_cells collect (each.H));
		P_moyenne <- mean(active_cells collect (each.P));
		temps_peche_moyen <- mean(pecheur collect (each.temps_peche));
		coeff_capture_jour <- biomasse_pechee_moyenne / ((H_moyenne + P_moyenne) * 0.5 * temps_peche_moyen); 
		
		create Mat_moyenne;
		save ("temps;nuit;ID_cell;grid_x;grid_y;pente_externe;statut;H_case_mensuel;CO_case_mensuel;P_case_mensuel;nb_pecheurs_case;biomasse_tot_pechee_herb_jour_case;biomasse_tot_pechee_herb_nuit_case;biomasse_tot_pechee_carniv_jour_case;biomasse_tot_pechee_carniv_nuit_case;biomasse_tot_pechee_herb_jour_case_r;biomasse_tot_pechee_herb_nuit_case_r;biomasse_tot_pechee_carniv_jour_case_r;biomasse_tot_pechee_carniv_nuit_case_r;nb_conflits_case_nuit;nb_conflits_case_jour_pecheur;nb_conflits_case_jour_tourisme;pecheurs_cell") to: output_csv_cell + ".txt" rewrite: false;
		save ("temps;ID_pecheur;ID_commune_pecheur;IDDist07_pecheur;radius;selectivite;dependance_ressource;poaching_probability;peche_nuit;temps_peche;proba_peche;nb_peche;nb_peche_pente_externe;nb_peche_lagon;nb_peche_AMP;nb_peche_hors_AMP;biomasse_tot_capturee_herb_lagon;biomasse_tot_capturee_herb_pente_externe;biomasse_tot_capturee_carniv_lagon;biomasse_tot_capturee_carniv_pente_externe;biomasse_tot_capturee_herb_lagon_r;biomasse_tot_capturee_herb_pente_externe_r;biomasse_tot_capturee_carniv_lagon_r;biomasse_tot_capturee_carniv_pente_externe_r;nb_conflits_pecheur_nuit;nb_conflits_pecheur_pecheur;nb_conflits_pecheur_tourisme;nb_quota_atteint_pecheur;cells_pecheur") to: output_csv_pecheur + ".txt" rewrite: false;
			
		
	} // Fermeture du init	
	
	//action qui permet de pre-calculer les distances entre ports possibles et cellules d'eau et qui sauvegarde le resultat dans un fichier 
	action pre_calcul_distances {
		list<cell> possible_cells;
		ask household_area {
			using topology(world) {possible_cells << active_cells with_min_of (each.location distance_to self);}
		}
		possible_cells <- remove_duplicates(possible_cells);
	    save "cell1,cell2,distance" to: chemin_fichier_distances;
	    int i <- 0;
	    int seuil <- 10;
	    ask possible_cells{
	    	list<cell> neigh;
	      	using topology(world) {
	      	  	neigh <- active_cells at_distance distance_bateau;
	      	}
	      	ask neigh {
	        	if (self = myself) {
		        	save ""+int(self)+","+int(myself)+",0" to: chemin_fichier_distances rewrite: false;
		        } else {
		        	path the_path <- active_cells path_between (self.location, myself.location);
					if (the_path != nil ) { // s'il n'y a pas de chemin jusqu'a la cellule ou qu'elle est trop loin, on l'enleve
						save ""+int(self)+","+int(myself)+","+round(the_path.shape.perimeter) to: chemin_fichier_distances rewrite: false;
					} else {
						save ""+int(self)+","+int(myself)+","+ "-1" to: chemin_fichier_distances rewrite: false;
					}	
				}
	        }
	        i <- i+1;
	        int avance <- round(i/length(possible_cells) * 100);
	        if (not mode_batch and (avance >= seuil)) {seuil <- seuil + 10; write "calcul des distances: " + avance +"%";}
	    }
	}
	
	//action qui permet de charger le fichier des distances
	action chargement_fichier_distances {
		matrix<int> ds <- matrix<int>(matrix(csv_file(chemin_fichier_distances, ",",int,true)));
		loop i from: 0 to: ds.rows - 1 {
			ask cell[ds[0,i]] {
				dists[ds[1,i]]<-float(ds[2,i]);
			}
			ask cell[ds[1,i]] {
				dists[ds[0,i]]<-float(ds[2,i]);
			}
		}
		if not mode_batch {write "chargement du fichier des distances fini ";}
	}
	
	//action qui permet de supprimer les cellule d'eau isolees (celles qui ne sont reliees qu'a peu d'autres cellules)
	action filtrage_cellules_eau_isolees {
		list<list<cell>> clusters <- list<list<cell>>(simple_clustering_by_distance(active_cells, 1));
        loop cluster over: clusters {
        	if (length(cluster) < seuil_filtrage_cellules_isolees) {
        		ask cluster {
    				color <- #orange;    		
	        		active_cells >> self;	
        		}
        	}
        }
	}
	
	
	// Différenciation jour - nuit
	reflex nuit when: nuit = true {      // Nuit	
	
		coeff_capture <- 2 * coeff_capture_jour;	
					
		ask active_cells parallel: true{
			tourism_criterion <- 0.0;     // Pas de tourisme la nuit		
		}
		
		ask pecheurs_de_jour {
			location <- my_house.location;  
			fishing_place <- nil;			
		}
		
		ask pecheurs_de_nuit {			
			if flip(proba_peche) {
				do choose_a_cell;  
       			do fishing;  
       		}	
       		else {
       			location <- my_house.location;  
				fishing_place <- nil;	
       		}				
		}		
	} 
	
	reflex jour when: nuit = false {      // Jour
	
		coeff_capture <- coeff_capture_jour;
				 
		ask active_cells parallel: true{
			tourism_criterion <- tourism_criterion_jour;
		}		
					
		ask pecheurs_de_nuit {
			location <- my_house.location;  
			fishing_place <- nil;					
		}
		
		ask pecheurs_de_jour {			
			if flip(proba_peche) {
				do choose_a_cell;  
       			do fishing;	       					
       		}	
       		else {
       			location <- my_house.location;  
				fishing_place <- nil;	
       		}	
		}
	} 
	
	reflex dynamic {		
		  	   
	   float max_biomass_poissons <- active_cells max_of (each.H + each.P);
	    	    
       ask active_cells parallel: true{
        	do dynamic;   
        	do compute_biomass_criterion(max_biomass_poissons);           	
       }
       
	   // Action de spillover quand on est proche de la carrying_capacity (90%)
	   ask active_cells where (each.H > (coeff_spill_over * each.carrying_capacity_H)) {
	    	do spillover_herbivores;	    	
	   }
	    
	    ask active_cells where (each.CO > (coeff_spill_over * each.carrying_capacity_CO)) {
	    	do spillover_corallivores;
	    }
	  
	    ask active_cells where (each.P > (coeff_spill_over * each.carrying_capacity_P)) {
	    	do spillover_predateurs;
	     }
	   
	    ask active_cells parallel: true{
	    	do biomass_update;	    	
	    } 
	}
	
	reflex pas_temps {
		ask active_cells {			
			nb_pecheur_present <- 0;
		}
	}
		
	reflex annees when: every (360.0 * 2) {		
		if temps > 0 {			
			ask active_cells {	
				H_case_annuel <- H_case_annuel / (360 * 2);
				P_case_annuel <- P_case_annuel / (360 * 2);
				biomasse_pechee_herb_jour_case_annuel_r <- biomasse_pechee_herb_jour_case_annuel_r / (360 * 2);
				biomasse_pechee_herb_nuit_case_annuel_r <- biomasse_pechee_herb_nuit_case_annuel_r / (360 * 2);			
				biomasse_pechee_carniv_jour_case_annuel_r <- biomasse_pechee_carniv_jour_case_annuel_r / (360 * 2);
				biomasse_pechee_carniv_nuit_case_annuel_r <- biomasse_pechee_carniv_nuit_case_annuel_r / (360 * 2);					
			}
			
			save active_cells with: [H_case_annuel :: "H", P_case_annuel :: "P", nb_pecheurs_case_annuel :: "nb_peche", biomasse_pechee_herb_jour_case_annuel :: "peche_h_j", biomasse_pechee_herb_nuit_case_annuel :: "peche_h_n", biomasse_pechee_carniv_jour_case_annuel :: "peche_c_j", biomasse_pechee_carniv_nuit_case_annuel :: "peche_c_n", biomasse_pechee_herb_jour_case_annuel_r :: "peche_hj_r", biomasse_pechee_herb_nuit_case_annuel_r :: "peche_hn_r", biomasse_pechee_carniv_jour_case_annuel_r :: "peche_cj_r", biomasse_pechee_carniv_nuit_case_annuel_r :: "peche_cn_r", nb_conflits_case_nuit_annuel :: "conflit_n", nb_conflits_case_jour_pecheur_annuel :: "conflit_p", nb_conflits_case_jour_tourisme_annuel :: "conflit_t"] to: output_carte_lagon + (temps / 360) + ".shp" type: "shp";			
			
			ask active_cells {
				H_case_annuel <- 0.0;
				P_case_annuel <- 0.0;
				biomasse_pechee_herb_jour_case_annuel <- 0.0;
				biomasse_pechee_herb_nuit_case_annuel <- 0.0;
				biomasse_pechee_carniv_jour_case_annuel <- 0.0;
				biomasse_pechee_carniv_nuit_case_annuel <- 0.0;
				biomasse_pechee_herb_jour_case_annuel_r <- 0.0;
				biomasse_pechee_herb_nuit_case_annuel_r <- 0.0;
				biomasse_pechee_carniv_jour_case_annuel_r <- 0.0;
				biomasse_pechee_carniv_nuit_case_annuel_r <- 0.0;
				nb_conflits_case_nuit_annuel <- 0;
				nb_conflits_case_jour_pecheur_annuel <- 0;
				nb_conflits_case_jour_tourisme_annuel <- 0;			
				nb_pecheurs_case_annuel <- 0;
			}
			
			ask commune {
				biomasse_pechee_herb_commune_annuel_r <- biomasse_pechee_herb_commune_annuel_r / (360 * 2);
				biomasse_pechee_carniv_commune_annuel_r <- biomasse_pechee_carniv_commune_annuel_r / (360 * 2);
			}
			
					
			save commune with: [ID_commune :: "commune", biomasse_pechee_herb_commune_annuel :: "peche_h", biomasse_pechee_carniv_commune_annuel :: "peche_c", biomasse_pechee_herb_commune_annuel_r :: "peche_h_r", biomasse_pechee_carniv_commune_annuel_r :: "peche_c_r", nb_conflits_nuit_commune_annuel :: "conflit_n", nb_conflits_jour_pecheur_commune_annuel :: "conflit_p", nb_conflits_jour_tourisme_commune_annuel :: "conflit_t", nb_quota_atteint_commune_annuel :: "quota"] to: output_carte_commune + (temps / 360) + ".shp" type: "shp";
		
			ask commune {
				biomasse_pechee_herb_commune_annuel <- 0.0;
				biomasse_pechee_carniv_commune_annuel <- 0.0;
				biomasse_pechee_herb_commune_annuel_r <- 0.0;
				biomasse_pechee_carniv_commune_annuel_r <- 0.0;
				nb_conflits_nuit_commune_annuel <- 0;
				nb_conflits_jour_pecheur_commune_annuel <- 0;
				nb_conflits_jour_tourisme_commune_annuel <- 0;
				nb_quota_atteint_commune_annuel <- 0;
			}
		}
	}
	
	// Fin de la pêche
			
	reflex save_results when: every (1 #cycle) {
		if (temps >= 0) {
			list<string> data;
		ask active_cells {
			//H_case_mensuel <- H_case_tot;
			//CO_case_mensuel <- CO_case_tot;
			//P_case_mensuel <- P_case_tot;
			//biomasse_tot_pechee_herb_jour_case_r <- biomasse_tot_pechee_herb_jour_case_r / (30 * 2);
			//biomasse_tot_pechee_herb_nuit_case_r <- biomasse_tot_pechee_herb_nuit_case_r / (30 * 2);
			//biomasse_tot_pechee_carniv_jour_case_r <- biomasse_tot_pechee_carniv_jour_case_r / (30 * 2);
			//biomasse_tot_pechee_carniv_nuit_case_r <- biomasse_tot_pechee_carniv_nuit_case_r / (30 * 2);
			
			data << ("" + temps + ";" + nuit + ";" + self + ";" + grid_x + ";" + grid_y + ";" + pente_externe + ";" + statut + ";" + H + ";" + CO + ";" + P + ";" + nb_pecheurs_case + ";" + biomasse_tot_pechee_herb_jour_case + ";" + biomasse_tot_pechee_herb_nuit_case + ";" + biomasse_tot_pechee_carniv_jour_case + ";" + biomasse_tot_pechee_carniv_nuit_case + ";" + biomasse_tot_pechee_herb_jour_case_r + ";" + biomasse_tot_pechee_herb_nuit_case_r + ";" + biomasse_tot_pechee_carniv_jour_case_r + ";" + biomasse_tot_pechee_carniv_nuit_case_r + ";" + nb_conflits_case_nuit + ";" + nb_conflits_case_jour_pecheur + ";" + nb_conflits_case_jour_tourisme + ";" + pecheurs_cell);
			//save ("" + temps + ";" + self + ";" + grid_x + ";" + grid_y + ";" + pente_externe + ";" + statut + ";" + H_case_tot + ";" + CO_case_tot + ";" + P_case_tot + ";" + nb_pecheurs_case + ";" + biomasse_tot_pechee_herb_jour_case + ";" + biomasse_tot_pechee_herb_nuit_case + ";" + biomasse_tot_pechee_carniv_jour_case + ";" + biomasse_tot_pechee_carniv_nuit_case + ";" + biomasse_tot_pechee_herb_jour_case_r + ";" + biomasse_tot_pechee_herb_nuit_case_r + ";" + biomasse_tot_pechee_carniv_jour_case_r + ";" + biomasse_tot_pechee_carniv_nuit_case_r + ";" + nb_conflits_case_nuit + ";" + nb_conflits_case_jour_pecheur + ";" + nb_conflits_case_jour_tourisme + ";" + pecheurs_cell) to: output_csv_cell + temps_annees + ".txt" rewrite: false;
	
			H_case_annuel <- H_case_annuel + H;
			P_case_annuel <- P_case_annuel + P;
			biomasse_pechee_herb_jour_case_annuel <- biomasse_pechee_herb_jour_case_annuel + biomasse_tot_pechee_herb_jour_case;
			biomasse_pechee_herb_nuit_case_annuel <- biomasse_pechee_herb_nuit_case_annuel + biomasse_tot_pechee_herb_nuit_case;			
			biomasse_pechee_carniv_jour_case_annuel <- biomasse_pechee_carniv_jour_case_annuel + biomasse_tot_pechee_carniv_jour_case;
			biomasse_pechee_carniv_nuit_case_annuel <- biomasse_pechee_carniv_nuit_case_annuel + biomasse_tot_pechee_carniv_nuit_case;
			biomasse_pechee_herb_jour_case_annuel_r <- biomasse_pechee_herb_jour_case_annuel_r + biomasse_tot_pechee_herb_jour_case_r;
			biomasse_pechee_herb_nuit_case_annuel_r <- biomasse_pechee_herb_nuit_case_annuel_r + biomasse_tot_pechee_herb_nuit_case_r;			
			biomasse_pechee_carniv_jour_case_annuel_r <- biomasse_pechee_carniv_jour_case_annuel_r + biomasse_tot_pechee_carniv_jour_case_r;
			biomasse_pechee_carniv_nuit_case_annuel_r <- biomasse_pechee_carniv_nuit_case_annuel_r + biomasse_tot_pechee_carniv_nuit_case_r;
			
			nb_conflits_case_nuit_annuel <- nb_conflits_case_nuit_annuel + nb_conflits_case_nuit;
			nb_conflits_case_jour_pecheur_annuel <- nb_conflits_case_jour_pecheur_annuel + nb_conflits_case_jour_pecheur;
			nb_conflits_case_jour_tourisme_annuel <- nb_conflits_case_jour_tourisme_annuel + nb_conflits_case_jour_tourisme;			
			nb_pecheurs_case_annuel <- nb_pecheurs_case_annuel + nb_pecheurs_case;
			
			nb_pecheurs_case <- 0; 
			biomasse_tot_pechee_herb_jour_case <- 0.0;
			biomasse_tot_pechee_herb_nuit_case <- 0.0;
			biomasse_tot_pechee_carniv_jour_case <- 0.0;
			biomasse_tot_pechee_carniv_nuit_case <- 0.0;
			biomasse_tot_pechee_herb_jour_case_r <- 0.0;
			biomasse_tot_pechee_herb_nuit_case_r <- 0.0;
			biomasse_tot_pechee_carniv_jour_case_r <- 0.0;
			biomasse_tot_pechee_carniv_nuit_case_r <- 0.0;
			nb_conflits_case_nuit <- 0;
			nb_conflits_case_jour_pecheur <- 0;
			nb_conflits_case_jour_tourisme <- 0;
			pecheurs_cell <- nil;
		}
		save data to: output_csv_cell + temps + ".txt" rewrite: false;
	
		ask pecheur {
			//biomasse_tot_capturee_herb_lagon_r <- biomasse_tot_capturee_herb_lagon_r / (30 * 2);
			//biomasse_tot_capturee_herb_pente_externe_r <- biomasse_tot_capturee_herb_pente_externe_r / (30 * 2);
			//biomasse_tot_capturee_carniv_lagon_r <- biomasse_tot_capturee_carniv_lagon_r / (30 * 2);
			//biomasse_tot_capturee_carniv_pente_externe_r <- biomasse_tot_capturee_carniv_pente_externe_r / (30 * 2);
		}
		ask commune {
			ask pecheur where (each.ID_commune_pecheur = ID_commune) {
				myself.biomasse_pechee_herb_commune_annuel <- myself.biomasse_pechee_herb_commune_annuel + biomasse_tot_capturee_herb_lagon + biomasse_tot_capturee_herb_pente_externe;
				myself.biomasse_pechee_carniv_commune_annuel <- myself.biomasse_pechee_carniv_commune_annuel + biomasse_tot_capturee_carniv_lagon + biomasse_tot_capturee_carniv_pente_externe;
				myself.biomasse_pechee_herb_commune_annuel_r <- myself.biomasse_pechee_herb_commune_annuel_r + biomasse_tot_capturee_herb_lagon_r + biomasse_tot_capturee_herb_pente_externe_r;
				myself.biomasse_pechee_carniv_commune_annuel_r <- myself.biomasse_pechee_carniv_commune_annuel_r + biomasse_tot_capturee_carniv_lagon_r + biomasse_tot_capturee_carniv_pente_externe_r;
				myself.nb_conflits_nuit_commune_annuel <- myself.nb_conflits_nuit_commune_annuel + nb_conflits_pecheur_nuit;
				myself.nb_conflits_jour_pecheur_commune_annuel <- myself.nb_conflits_jour_pecheur_commune_annuel + nb_conflits_pecheur_pecheur;
				myself.nb_conflits_jour_tourisme_commune_annuel <- myself.nb_conflits_jour_tourisme_commune_annuel + nb_conflits_pecheur_tourisme;
				myself.nb_quota_atteint_commune_annuel <- myself.nb_quota_atteint_commune_annuel + nb_quota_atteint_pecheur;
			}
			
		}
		data <- [];
		
		ask pecheur {
			
			//save [temps, self, IDDist07_pecheur, selectivite, poaching_probability, peche_nuit, temps_peche, proba_peche, nb_peche, nb_peche_pente_externe, nb_peche_lagon, nb_peche_AMP, nb_peche_hors_AMP, biomasse_tot_capturee_herb_lagon, biomasse_tot_capturee_herb_pente_externe, biomasse_tot_capturee_carniv_lagon, biomasse_tot_capturee_carniv_pente_externe, nb_conflits_pecheur_nuit, nb_conflits_pecheur_pecheur, nb_conflits_pecheur_tourisme, nb_quota_atteint_pecheur, cells_pecheur] to: output_csv_pecheur + temps_annees + ".csv" type:"csv" rewrite: false;
			
			//save ("" + temps + ";" + self + ";" + ID_commune_pecheur + ";" + IDDist07_pecheur + ";" + fishing_radius + ";" + selectivite + ";" + dependance_ressource_pecheur + ";" + poaching_probability + ";" + peche_nuit + ";" + temps_peche + ";" + proba_peche + " ; " + nb_peche + ";" + nb_peche_pente_externe + ";" + nb_peche_lagon + ";" + nb_peche_AMP + ";" + nb_peche_hors_AMP + ";" + biomasse_tot_capturee_herb_lagon + ";" + biomasse_tot_capturee_herb_pente_externe + ";" + biomasse_tot_capturee_carniv_lagon + " ; " + biomasse_tot_capturee_carniv_pente_externe + " ; " + biomasse_tot_capturee_herb_lagon_r + ";" + biomasse_tot_capturee_herb_pente_externe_r + ";" + biomasse_tot_capturee_carniv_lagon_r + " ; " + biomasse_tot_capturee_carniv_pente_externe_r + " ; " + nb_conflits_pecheur_nuit + ";" + nb_conflits_pecheur_pecheur + ";" + nb_conflits_pecheur_tourisme + ";" + nb_quota_atteint_pecheur + ";" + cells_pecheur) to: output_csv_pecheur + temps_annees + ".csv" type:"text" rewrite: false ;
			data << ("" + temps + ";" + self + ";" + ID_commune_pecheur + ";" + IDDist07_pecheur + ";" + fishing_radius + ";" + selectivite + ";" + dependance_ressource_pecheur + ";" + poaching_probability + ";" + peche_nuit + ";" + temps_peche + ";" + proba_peche + " ; " + nb_peche + ";" + nb_peche_pente_externe + ";" + nb_peche_lagon + ";" + nb_peche_AMP + ";" + nb_peche_hors_AMP + ";" + biomasse_tot_capturee_herb_lagon + ";" + biomasse_tot_capturee_herb_pente_externe + ";" + biomasse_tot_capturee_carniv_lagon + " ; " + biomasse_tot_capturee_carniv_pente_externe + " ; " + biomasse_tot_capturee_herb_lagon_r + ";" + biomasse_tot_capturee_herb_pente_externe_r + ";" + biomasse_tot_capturee_carniv_lagon_r + " ; " + biomasse_tot_capturee_carniv_pente_externe_r + " ; " + nb_conflits_pecheur_nuit + ";" + nb_conflits_pecheur_pecheur + ";" + nb_conflits_pecheur_tourisme + ";" + nb_quota_atteint_pecheur + ";" + cells_pecheur);
							
			biomasse_tot_capturee_herb_lagon <- 0.0;
			biomasse_tot_capturee_herb_pente_externe <- 0.0;
			biomasse_tot_capturee_carniv_lagon <- 0.0;
			biomasse_tot_capturee_carniv_pente_externe <- 0.0;
			biomasse_tot_capturee_herb_lagon_r <- 0.0;
			biomasse_tot_capturee_herb_pente_externe_r <- 0.0;
			biomasse_tot_capturee_carniv_lagon_r <- 0.0;
			biomasse_tot_capturee_carniv_pente_externe_r <- 0.0;
			nb_conflits_pecheur_nuit <- 0;
			nb_conflits_pecheur_pecheur <- 0;
			nb_conflits_pecheur_tourisme <- 0;
			nb_quota_atteint_pecheur <- 0;	
			nb_peche <- 0;
			nb_peche_AMP <- 0;
			nb_peche_hors_AMP <- 0;
			nb_peche_pente_externe <- 0;
			nb_peche_lagon <- 0;
			cells_pecheur <- nil;
		}
		save data to: output_csv_pecheur + temps + ".txt" rewrite: false;
		}	
	}
	
	reflex fin_journee {
				
		H_moyenne <- mean(active_cells collect (each.H));
		CO_moyenne <- mean(active_cells collect (each.CO));
		P_moyenne <- mean(active_cells collect (each.P));
		H_moyenne_lagon <- active_cells where (each.pente_externe = false) mean_of (each.H);
		CO_moyenne_lagon <- active_cells where (each.pente_externe = false) mean_of (each.CO);
		P_moyenne_lagon <- active_cells where (each.pente_externe = false) mean_of (each.P);
		H_moyenne_pente_externe <- active_cells where (each.pente_externe = true) mean_of (each.H);
		CO_moyenne_pente_externe <- active_cells where (each.pente_externe = true) mean_of (each.CO);  
		P_moyenne_pente_externe <- active_cells where (each.pente_externe = true) mean_of (each.P); 
		C_moyenne_pente_externe <- active_cells where (each.pente_externe = true) mean_of (each.C); 
		C_moyenne_lagon <- active_cells where (each.pente_externe = false) mean_of (each.C);
		T_moyenne_pente_externe <- active_cells where (each.pente_externe = true) mean_of (each.T); 
		T_moyenne_lagon <- active_cells where (each.pente_externe = false) mean_of (each.T);
		
		H_moyenne_AMP <- mean(active_cells where (each.statut = 1) collect (each.H));
		CO_moyenne_AMP <- mean(active_cells where (each.statut = 1) collect (each.CO));
		P_moyenne_AMP <- mean(active_cells where (each.statut = 1) collect (each.P));
		H_moyenne_not_AMP <- mean(active_cells where (each.statut = 0) collect (each.H));
		CO_moyenne_not_AMP <- mean(active_cells where (each.statut = 0) collect (each.CO));
		P_moyenne_not_AMP <- mean(active_cells where (each.statut = 0) collect (each.P));
		
		biomasse_pechee_herb_pente_externe <- active_cells where (each.pente_externe = true) sum_of (each.biomasse_pechee_herb_case); 
		biomasse_pechee_herb_lagon <- active_cells where (each.pente_externe = false) sum_of (each.biomasse_pechee_herb_case);
		biomasse_pechee_carniv_pente_externe <- active_cells where (each.pente_externe = true) sum_of (each.biomasse_pechee_carniv_case); 
		biomasse_pechee_carniv_lagon <- active_cells where (each.pente_externe = false) sum_of (each.biomasse_pechee_carniv_case); 
		biomasse_pechee_herb_pente_externe_r <- active_cells where (each.pente_externe = true) mean_of (each.biomasse_pechee_herb_case_r); 
		biomasse_pechee_herb_lagon_r <- active_cells where (each.pente_externe = false) mean_of (each.biomasse_pechee_herb_case_r);
		biomasse_pechee_carniv_pente_externe_r <- active_cells where (each.pente_externe = true) mean_of (each.biomasse_pechee_carniv_case_r); 
		biomasse_pechee_carniv_lagon_r <- active_cells where (each.pente_externe = false) mean_of (each.biomasse_pechee_carniv_case_r); 
		
		conflits_nuit <- active_cells sum_of (each.conflits_case_nuit);		
		conflits_jour_tourisme <- active_cells sum_of (each.conflits_case_jour_tourisme);
		conflits_jour_pecheur <- active_cells sum_of (each.conflits_case_jour_pecheur);
		
		quota_atteint <- pecheur sum_of (each.quota_atteint_pecheur);
		
		ask Mat_moyenne {   // Matrice pour csv
			temps_M <- temps;
			nuit_M <- nuit;
			H_moyenne_pente_externe_M <- H_moyenne_pente_externe;
			H_moyenne_lagon_M <- H_moyenne_lagon;			
			CO_moyenne_M <- CO_moyenne;
			P_moyenne_pente_externe_M <- P_moyenne_pente_externe;
			P_moyenne_lagon_M <- P_moyenne_lagon;
			C_moyenne_pente_externe_M <- C_moyenne_pente_externe; 
			C_moyenne_lagon_M <- C_moyenne_lagon;
			T_moyenne_pente_externe_M <- T_moyenne_pente_externe; 
			T_moyenne_lagon_M <- T_moyenne_lagon;
			biomasse_pechee_herb_pente_externe_M <- biomasse_pechee_herb_pente_externe;
			biomasse_pechee_herb_lagon_M <- biomasse_pechee_herb_lagon;
			biomasse_pechee_carniv_pente_externe_M <- biomasse_pechee_carniv_pente_externe;
			biomasse_pechee_carniv_lagon_M <- biomasse_pechee_carniv_lagon;
			biomasse_pechee_herb_pente_externe_M_r <- biomasse_pechee_herb_pente_externe_r;
			biomasse_pechee_herb_lagon_M_r <- biomasse_pechee_herb_lagon_r;
			biomasse_pechee_carniv_pente_externe_M_r <- biomasse_pechee_carniv_pente_externe_r;
			biomasse_pechee_carniv_lagon_M_r <- biomasse_pechee_carniv_lagon_r;
			
			conflits_nuit_M <- conflits_nuit;			
			conflits_jour_pecheur_M <- conflits_jour_pecheur;	
			conflits_jour_tourisme_M <- conflits_jour_tourisme;
			quota_atteint_M <- quota_atteint;
			
			save Mat_moyenne to: output_csv_ecologique type: "csv" rewrite: false;
		}		
						
			
		temps <- temps + 0.5;
		
		ask pecheur {
			biomasse_capturee_herb <- 0.0;
			biomasse_capturee_carniv <- 0.0;
			biomasse_capturee_herb_r <- 0.0;
			biomasse_capturee_carniv_r <- 0.0;
			quota_atteint_pecheur <- 0;			
		}		
		ask active_cells {
			H_pechee <- 0.0;
			P_pechee <- 0.0;
			pecheur_present <- false;
			pecheur_present_criterion <- 0.4; // Pas de pêche présent -> bonne préférence pour la case
			conflits_case_nuit <- 0;
			conflits_case_jour_pecheur <- 0;
			conflits_case_jour_tourisme <- 0;	
			biomasse_pechee_herb_case <- 0.0;
			biomasse_pechee_carniv_case <- 0.0;	
			biomasse_pechee_herb_case_r <- 0.0;
			biomasse_pechee_carniv_case_r <- 0.0;		
		}
		
		if nuit {
			nuit <- false;
		}
		else {
			nuit <- true;
		}
	}	
	
		
	reflex when: mode_batch {
		save [cycle,H_moyenne_lagon,CO_moyenne_lagon,P_moyenne_lagon,H_moyenne_pente_externe,CO_moyenne_pente_externe,P_moyenne_pente_externe,			
			C_moyenne_lagon,T_moyenne_lagon,
			C_moyenne_pente_externe,T_moyenne_pente_externe] type: "csv" to: resultats_batch rewrite: false;
	}
	// Fin de la simulation à t = temps_simulation ; sorties d'outputs	
	reflex end_simulation when: temps = temps_simulation {
		do pause;
		
	}
	
}   // Fermeture du global

species Mat_moyenne {	
	float temps_M;
	bool nuit_M;
	float H_moyenne_pente_externe_M;
	float H_moyenne_lagon_M;
	float CO_moyenne_M;
	float P_moyenne_pente_externe_M;
	float P_moyenne_lagon_M;
	float C_moyenne_pente_externe_M; 
	float C_moyenne_lagon_M;
	float T_moyenne_pente_externe_M; 
	float T_moyenne_lagon_M;
	float biomasse_pechee_herb_pente_externe_M;
	float biomasse_pechee_herb_lagon_M;
	float biomasse_pechee_carniv_pente_externe_M;
	float biomasse_pechee_carniv_lagon_M;
	float biomasse_pechee_herb_pente_externe_M_r;
	float biomasse_pechee_herb_lagon_M_r;
	float biomasse_pechee_carniv_pente_externe_M_r;
	float biomasse_pechee_carniv_lagon_M_r;
	float nb_conflits_M;		
	int conflits_nuit_M;
	int conflits_jour_pecheur_M;
	int	conflits_jour_tourisme_M;
	int quota_atteint_M;
}

species Mat_cell {	
	float temps_t;
	string periode;
	float x_location;
	float y_location;
	string pente_externe_cell;
	string AMP;
	float H_cell;
	float CO_cell;
	float P_cell;
	float C_cell;
	float T_cell;
	float ID_pecheur_cell;
	float biom_pechee_herb_cell;
	float biom_pechee_carniv_cell;
	int conflit_nuit_cell;
	int conflit_pecheur_cell;
	int conflit_tourisme_cell;
	int quota_atteint_cell;
}


species pecheur {
	
    bool en_train_pecher <- false;
	household_area my_house <- nil;
	list<cell> fishing_area;
	float fishing_radius;
	cell fishing_place_tmp;
	cell fishing_place;
	map<cell,float> distances;	
	list<cell> fishing_area_tot;
		
	float poaching_probability;
	float dependance_ressource_pecheur;
	float selectivite;
	bool peche_nuit;
	float proba_peche;
	bool bateau_pecheur <- false;
	bool pirogue_pecheur <- false;	
	float temps_peche; // en heure
	
	float biomasse_capturee_herb;
	float biomasse_tot_capturee_herb_lagon;
	float biomasse_tot_capturee_herb_pente_externe;
	float biomasse_capturee_carniv;	
	float biomasse_tot_capturee_carniv_lagon;
	float biomasse_tot_capturee_carniv_pente_externe;
	float biomasse_capturee_herb_r;
	float biomasse_tot_capturee_herb_lagon_r;
	float biomasse_tot_capturee_herb_pente_externe_r;
	float biomasse_capturee_carniv_r;	
	float biomasse_tot_capturee_carniv_lagon_r;
	float biomasse_tot_capturee_carniv_pente_externe_r;
	int nb_conflits_pecheur_nuit;
	int nb_conflits_pecheur_pecheur;
	int nb_conflits_pecheur_tourisme;
	int quota_atteint_pecheur;	
	int nb_quota_atteint_pecheur;
	float nb_quota_atteint_pecheur_annuels;
	float nb_quota_pecheur_annuel;
	int IDDist07_pecheur;	
	int ID_LOG_pecheur;
	int ID_commune_pecheur;
	
	int nb_peche;	
	int nb_peche_AMP;
	int nb_peche_hors_AMP;
	int nb_peche_lagon;
	int nb_peche_pente_externe;
	list<cell> cells_pecheur;	
	
	cell cell_pecheur; //cellule d'eau la plus proche de moi
	float dist_to_cell_pecheur;			
	  			
	action init_pecheur {
		
		ask District where (each.IDDist07 = IDDist07_pecheur) {
			myself.poaching_probability <- sensitivity_district / (1 + surveillance);			
	  		myself.dependance_ressource_pecheur <- sensitivity_district;
	  		myself.ID_commune_pecheur <- ID_commune_district;
			my_households <- household_area overlapping self;
			myself.my_house <- one_of(my_households);
		}
		if (my_house != nil) {
	       	location <- my_house.location;         
	    }
	    else {
	     	do die;
	    }     
		
		if flip(0.7) {
			if interdiction_peche_nocturne = false {
				peche_nuit <- true;     // Pêche uniquement de nuit
				poaching_probability <- 1.0 / (1 + surveillance); 							
			}
			else {
				peche_nuit <- false;																	
			}
		}
		else {					
			peche_nuit <- false;  // Pêche uniquement de jour								
		}	
						  				
	  	selectivite <- rnd(100)/100;
	  	if selectivite = 1 {
	  		selectivite <- 0.99;
	  	}	
		
		using topology(world) {
			cell_pecheur <- active_cells with_min_of (each.location distance_to self);
			dist_to_cell_pecheur <- cell_pecheur distance_to self;
		}
		location <- cell_pecheur.location; //on deplace le pecheur sur sa cellule ou il y a son bateau
				                        		
	  using topology(world) {
			fishing_area_tot <- active_cells at_distance fishing_radius;
		}
		list<cell> cells_to_remove <- []; //liste des cellules inaccessibles
		loop fa over: fishing_area_tot
		{
			//si on est deja sur la cellule, distance de 0 et donc critere distance a 1
			if (fa = cell_pecheur) {
				distances[fa] <- 1;
			//sinon on cherche dans les distances deja enregistrees
			} else if (int(fa) in cell_pecheur.dists.keys) {
				float dist <- cell_pecheur.dists[int(fa)];
				if (dist < 0 or (dist > fishing_radius)) {
					cells_to_remove << fa;
				} else {	
					distances[fa] <- 1 - (ln((dist / fishing_radius) + 1) / ln(2));
				}
			//si la distance n'a pas encore ete calculee, on la calcule		
			} else {
				path the_path <- fishing_area_tot path_between (cell_pecheur.location, fa.location);
				if (the_path = nil or (the_path.shape.perimeter > fishing_radius))
				{ // s'il n'y a pas de chemin jusqu'a la cellule ou qu'elle est trop loin, on l'enleve
					cells_to_remove << fa;
					cell_pecheur.dists[int(fa)] <- -1.0;
				} else {
					float dist <- the_path.shape.perimeter;
					cell_pecheur.dists[int(fa)] <- dist;
					distances[fa] <- 1 - (ln((dist / fishing_radius) + 1) / ln(2));
				}	
			}
		}
		fishing_area_tot <- fishing_area_tot - cells_to_remove;
	    
	    if (flip(poaching_probability)){
	    	fishing_area <- fishing_area_tot ; // En cas de braconnage, pêche aussi dans les AMPs
	    }
	    else {
	    	fishing_area <- fishing_area_tot where not (each.statut = 1);  // Pas de pêche dans les AMPs si pas de braconnage
	    } 		    		
		location <- my_house.location;     //on replace le pecheur dans sa maison
	}
	
	action choose_a_cell {		
		
		// Le choix de la cellule se fait en fonction du critère de préférence, du tourisme, de la distance à parcourir et de la biomasse de poissons des cellules
		list<float> vals <- fishing_area collect (each.biomass_criterion + each.pref_criterion + each.tourism_criterion + distances[each] + each.pecheur_present_criterion); 
		fishing_place <- fishing_area[rnd_choice(vals)];
		
		if ((fishing_place != nil) and fishing_place.pecheur_present = true) {	
			fishing_place_tmp <- fishing_place;
			
			ask fishing_place {	
				if nuit {
					conflits_case_nuit <- conflits_case_nuit + 1;
					nb_conflits_case_nuit <- nb_conflits_case_nuit + 1;					
	    			myself.nb_conflits_pecheur_nuit <- myself.nb_conflits_pecheur_nuit + 1; 	    						
				}	
				else {
					conflits_case_jour_pecheur <- conflits_case_jour_pecheur + 1;	
					nb_conflits_case_jour_pecheur <- nb_conflits_case_jour_pecheur + 1;
	    			myself.nb_conflits_pecheur_pecheur <- myself.nb_conflits_pecheur_pecheur + 1;	
				}
						
			}									
									
			fishing_place <- one_of(fishing_place.neighbors where (each.pecheur_present = false) inter fishing_area);	
			if fishing_place = nil {
				fishing_place <- one_of(fishing_place_tmp.neighbors inter fishing_area);								
			}				
		}		
		
		if (fishing_place != nil) {
			location <- fishing_place.location;	
		
			ask fishing_place {
				
				add myself to: pecheurs_cell;
				add self to: myself.cells_pecheur;
				
				nb_pecheurs_case <- nb_pecheurs_case + 1;
				myself.nb_peche <- myself.nb_peche + 1;
				if statut = 1 {
					myself.nb_peche_AMP <- myself.nb_peche_AMP + 1;
				}
				else {
					myself.nb_peche_hors_AMP <- myself.nb_peche_hors_AMP + 1;
				}
				
				if (pente_externe) {
					myself.nb_peche_pente_externe <- myself.nb_peche_pente_externe + 1;					
				}
				else {
					myself.nb_peche_lagon <- myself.nb_peche_lagon + 1;		
				}
				
				float p <- rnd (0,100) / 100;
				if (p < level_tourism) and (nuit = false) {
					conflits_case_jour_tourisme <- conflits_case_jour_tourisme + 1;
					nb_conflits_case_jour_tourisme <- nb_conflits_case_jour_tourisme + 1;
					myself.nb_conflits_pecheur_tourisme <- myself.nb_conflits_pecheur_tourisme + 1;	    			
				}
				pecheur_present <- true;
				nb_pecheur_present <- nb_pecheur_present + 1;
				pecheur_present_criterion <- 0.2;	// pêcheur présent -> moins de préférence pour la case	
			}
			
		}		
	}	
	
	action fishing {
					
		if (fishing_place != nil) {
		// Biomasse de poissons pêchée par un pêcheur à une sortie		
		ask fishing_place {			
			
			float fishing_mass <-  (H + P) * myself.temps_peche * (1 - myself.selectivite) * coeff_capture;
					
			if (H != 0.0 or P != 0) {   // Pêche des herbivores
				if (fishing_mass >= quota) {
					myself.biomasse_capturee_herb <- min([(fishing_mass * H) / (H + P), (quota * H) / (H + P), H]);	
					myself.biomasse_capturee_carniv <- min([(fishing_mass * P) / (H + P), (quota * P) / (H + P), P]);
	    			myself.quota_atteint_pecheur <- myself.quota_atteint_pecheur + 1;
	    			myself.nb_quota_atteint_pecheur <- myself.nb_quota_atteint_pecheur + 1; 
				}
				else {	
					myself.biomasse_capturee_herb <- min([(fishing_mass * H) / (H + P), H]);				
					myself.biomasse_capturee_carniv <- min([(fishing_mass * P) / (H + P), P]);
				}
				
				myself.biomasse_capturee_herb_r <- myself.biomasse_capturee_herb / H;
				myself.biomasse_capturee_carniv_r <- myself.biomasse_capturee_carniv / P;   					   		
			} 
			
			H_pechee <- myself.biomasse_capturee_herb / 8;
			P_pechee <- myself.biomasse_capturee_carniv / 5;		

			biomasse_pechee_herb_case <- biomasse_pechee_herb_case + myself.biomasse_capturee_herb;
			biomasse_pechee_herb_case_r <- biomasse_pechee_herb_case_r + myself.biomasse_capturee_herb_r;
			biomasse_pechee_carniv_case <- biomasse_pechee_carniv_case + myself.biomasse_capturee_carniv;	
			biomasse_pechee_carniv_case_r <- biomasse_pechee_carniv_case_r + myself.biomasse_capturee_carniv_r;	
		   	
		   	if nuit {
		   		biomasse_tot_pechee_herb_nuit_case <- biomasse_tot_pechee_herb_nuit_case + myself.biomasse_capturee_herb; 
		   		biomasse_tot_pechee_herb_nuit_case_r <- biomasse_tot_pechee_herb_nuit_case_r + myself.biomasse_capturee_herb_r; 
				biomasse_tot_pechee_carniv_nuit_case <- biomasse_tot_pechee_carniv_nuit_case + myself.biomasse_capturee_carniv;
				biomasse_tot_pechee_carniv_nuit_case_r <- biomasse_tot_pechee_carniv_nuit_case_r + myself.biomasse_capturee_carniv_r;  			
			}
			else {
				biomasse_tot_pechee_herb_jour_case <- biomasse_tot_pechee_herb_jour_case + myself.biomasse_capturee_herb; 
				biomasse_tot_pechee_herb_jour_case_r <- biomasse_tot_pechee_herb_jour_case_r + myself.biomasse_capturee_herb_r; 
				biomasse_tot_pechee_carniv_jour_case <- biomasse_tot_pechee_carniv_jour_case + myself.biomasse_capturee_carniv; 
				biomasse_tot_pechee_carniv_jour_case_r <- biomasse_tot_pechee_carniv_jour_case_r + myself.biomasse_capturee_carniv_r; 
			}
			if pente_externe {			
				myself.biomasse_tot_capturee_herb_pente_externe <- myself.biomasse_tot_capturee_herb_pente_externe + myself.biomasse_capturee_herb;
				myself.biomasse_tot_capturee_herb_pente_externe_r <- myself.biomasse_tot_capturee_herb_pente_externe_r + myself.biomasse_capturee_herb_r;
				myself.biomasse_tot_capturee_carniv_pente_externe <- myself.biomasse_tot_capturee_carniv_pente_externe + myself.biomasse_capturee_carniv;
				myself.biomasse_tot_capturee_carniv_pente_externe_r <- myself.biomasse_tot_capturee_carniv_pente_externe_r + myself.biomasse_capturee_carniv_r;
			}
			else {
				myself.biomasse_tot_capturee_herb_lagon <- myself.biomasse_tot_capturee_herb_lagon + myself.biomasse_capturee_herb;
				myself.biomasse_tot_capturee_herb_lagon_r <- myself.biomasse_tot_capturee_herb_lagon_r + myself.biomasse_capturee_herb_r;
				myself.biomasse_tot_capturee_carniv_lagon <- myself.biomasse_tot_capturee_carniv_lagon + myself.biomasse_capturee_carniv;
				myself.biomasse_tot_capturee_carniv_lagon_r <- myself.biomasse_tot_capturee_carniv_lagon_r + myself.biomasse_capturee_carniv_r;
			}
		   				
		}
			
		}
	 }
		
	aspect base {
		draw circle (20) color: #magenta;
	}
	
} // Fermeture de pêcheur


// Grille d'unite spatiale
grid cell file: coral_file frequency: 0 neighbors: 8 {
	
	bool pente_externe <- false;
	bool clipped <- false;
	bool is_active <- false;
	bool active <- true;
	int statut;
	bool phase_perturbation <- false;
	
	// Variables du modèle écologique
	float C; // Biomasse de corail
	float T;   // Biomasse d'algues
	float M;
	float H;    // Biomasse d'herbivores
	float H_tmp;
	float CO;    // Biomasse de corallivores
	float CO_tmp;
	float P;      // Biomasse de prédateurs
	float P_tmp;
	float C_init;
	float T_init;
	float H_init;
	float P_init;
	float CO_init;
	float carrying_capacity_H;				
	float carrying_capacity_CO;		
	float carrying_capacity_P;
	float carrying_capacity_C_T;	
	float H_min;
	float CO_min;
	float P_min;
	float beta_C_CO_2;
	float beta_C_T_2;
	float beta_C_A_2;
	float beta_T_C_2;
	float beta_T_H_2;	
	float beta_H_P_2;
	float delta_H_T_2;
	float delta_H_A_2;
	float beta_CO_P_2;
	float delta_CO_C_2;
	float delta_P_H_CO_2;	
	float beta_H_PE;
	float beta_P_PE;
	float alpha_C;
	float alpha_T;	
	float H_pechee;
	float P_pechee;
	float biomasse_pechee_herb_case;
	float biomasse_pechee_carniv_case;
	float biomasse_tot_pechee_herb_jour_case;
	float biomasse_tot_pechee_herb_nuit_case;
	float biomasse_pechee_herb_jour_case_annuel;
	float biomasse_pechee_herb_nuit_case_annuel;
	float biomasse_tot_pechee_carniv_jour_case;
	float biomasse_tot_pechee_carniv_nuit_case;
	float biomasse_pechee_carniv_jour_case_annuel;
	float biomasse_pechee_carniv_nuit_case_annuel;
	float biomasse_pechee_herb_case_r;
	float biomasse_pechee_carniv_case_r;
	float biomasse_tot_pechee_herb_jour_case_r;
	float biomasse_tot_pechee_herb_nuit_case_r;
	float biomasse_pechee_herb_jour_case_annuel_r;
	float biomasse_pechee_herb_nuit_case_annuel_r;
	float biomasse_tot_pechee_carniv_jour_case_r;
	float biomasse_tot_pechee_carniv_nuit_case_r;
	float biomasse_pechee_carniv_jour_case_annuel_r;
	float biomasse_pechee_carniv_nuit_case_annuel_r;
	int conflits_case_nuit;
	int nb_conflits_case_nuit;
	int nb_conflits_case_nuit_annuel;
	int conflits_case_jour_tourisme;
	int nb_conflits_case_jour_tourisme_annuel;
	int nb_conflits_case_jour_tourisme;
	int conflits_case_jour_pecheur;
	int nb_conflits_case_jour_pecheur;
	int nb_conflits_case_jour_pecheur_annuel;
	float H_case_annuel;
	float CO_case_annuel;
	float P_case_annuel;
	int nb_pecheurs_case;
	int nb_pecheurs_case_annuel;
	
	int COTS <- 0;
	bool phase_recuperation <- false;
	bool eq_atteint <- false; 		
	bool perturbation_passee <- false;
	float recuperation;
	float temps_recuperation_corail_case;
		
	bool pecheur_present <- false;	
	float level_tourism;
	float tourism_criterion;
	float tourism_criterion_jour;
	float preference;
	float pref_criterion;
	float biomass_criterion;
	float pecheur_present_criterion;	
	int nb_pecheur_present;
		
	map<int,float> dists;
	list<pecheur> pecheurs_cell;
					
	action compute_biomass_criterion(float max_biomass_poissons) {		
		if (H + P) > 0 {
			biomass_criterion <- ln ((H + P) / max_biomass_poissons + 1) / ln(2); // Critère de biomasse lissé
		}
		else {
			biomass_criterion <- 0.0;
		}	
	}	
	
		
	// Modèle écologique de Lotka Voltera
	
	action dynamic {    
				        
        // Réinitialisation des paramètres de Lotka Volterra quand on atteint de nouveau l'état initial / d'équilibre
        if (perturbation_passee = true and eq_atteint = false) {      							
			eq_atteint <- true;
			recuperation <- temps;
			phase_perturbation <- false;      		
     	 } 
     	   
            	       
        C <- C + step * C * (alpha_C - coeff_destruction_corail_COTS * COTS - beta_C_CO_2 * CO - beta_C_T_2 * T);        	 		
     	T <- T + step * T * (alpha_T - beta_T_H_2 * H - beta_T_C_2 * C); 
     	
     	H_tmp <- H + step * H * (alpha_H - gamma_H - beta_H_P_2 * P - beta_H_PE * nb_pecheur_present + delta_H_T_2 * T) - H_pechee; 		  		  		
    	CO_tmp <- CO + step * CO * (alpha_CO - gamma_CO - beta_CO_P_2 * P + delta_CO_C_2 * C);    	
     	P_tmp <- P + step * P * (alpha_P - gamma_P - beta_P_PE * nb_pecheur_present + delta_P_H_CO_2 * (H + CO)) - P_pechee;      	
      	
      	if pente_externe {
        if (temps > debut_perturbation and perturbation = true and perturbation_passee = false) {   // perturbation
       		phase_perturbation <- true;
       		if C <= C_min {
       			phase_recuperation <- true;
       			COTS <- 0;   
       			alpha_C <- alpha_C * 200;       			
       		}      	
      	
      		if phase_recuperation = false { // Phase 1      			
      			COTS <- 1;       		
      		}  
      		else {  // Phase 2      		
      			if (C >= C_init) {
      				perturbation_passee <- true;
       			}
       			//write("phase recuperation");
       			
      		}
      	} 
      }
      	
      	
      	if (C >= C_init){      	
      		C <- C_init;
     	} 
     	if (T <= T_init) {
     		T <- T_init;
     	}  
     	
     	if C <= C_min {
			C <- C_min;				
		}	
			
		if T <= T_min {
			T <- T_min;				
		}	
     	
     	if (C + T >= carrying_capacity_C_T) {
     		T <- carrying_capacity_C_T - C;     		
     	}  	
      
     	if (H_tmp <= H_min) {
     		H_tmp <- H_min;
     	}
     	if (CO_tmp <= CO_min) {
     	 	CO_tmp <- CO_min;
     	}
     	if (P_tmp <= P_min) { 		
     		P_tmp <- P_min;
     	}
     		
    }
    
  	
	action spillover_herbivores {    // Spillover : déplacement de 10% des poissons présent dans la cellule 
		
	 	float quantity_herbivores <- (H - (carrying_capacity_H * coeff_spill_over));
		ask neighbors {			
	    	H_tmp <- H_tmp + quantity_herbivores / length(neighbors);	    	
	    }	    
	   		H_tmp <- H_tmp - quantity_herbivores;
	}
	
	action spillover_corallivores {
		float quantity_corallivores <- (CO -(carrying_capacity_CO * coeff_spill_over)); 
		ask neighbors {			
	    	CO_tmp <- CO_tmp + quantity_corallivores / length(neighbors);	    	
	    }	    
	   		CO_tmp <- CO_tmp - quantity_corallivores;	   	
	}
	
	action spillover_predateurs {
		
		float quantity_predateurs <- (P -(carrying_capacity_P * coeff_spill_over)); 
		ask neighbors {			
	    	P_tmp <- P_tmp + quantity_predateurs/length(neighbors);	    	
	    }	    
	   		P_tmp <- P_tmp - quantity_predateurs;
	}	
	
	action biomass_update {  // mise à jour à cause du spill over
		H <- H_tmp; 
		CO <- CO_tmp;
		P <- P_tmp;
	}
	
	aspect default {
		if (is_active) {
			draw shape color: rgb(C * 255, T * 255, 0);	
		}
		if (pente_externe) {
			draw pyramid(50) color: #yellow border: #black;
		}
	}
}


species District {
	
	list<household_area> my_households;
	list<pecheur> my_fishers;
	int IDDist07;
	int ID_commune_district;
	string bord_mer;
	float sensitivity_district;
	int nb_pecheurs_annexes_district;
	int nb_pecheurs_pro_district;
	int nb_individus_district;	
	float biomasse_pechee_herb_district_annuel;
	float biomasse_pechee_carniv_district_annuel;
	int nb_conflits_nuit_district_annuel;
	int nb_conflits_jour_pecheur_district_annuel;
	int nb_conflits_jour_tourisme_district_annuel;
	int nb_quota_atteint_district_annuel;
			
	aspect basic{
		color <- #lightgrey;
		draw shape color: color border: #grey;	
	}
}

species amp {
	int statut_AMP <- 1;
	aspect default {
		draw shape color:statut_AMP = 1 ? #red : #white border: #black;
	}
}

species household {
	int ID_District;
	int ID_LOG;
    int bateaux;
    int pirogues;    
    int pecheurs_pro;
    int pecheurs_annexes;    
}

species household_area {
    aspect base {
        draw triangle(100) color: #green border: #black;
    }
}

species commune {
	int ID_commune;
	int nb_pecheurs_annexes_commune;
	int nb_pecheurs_pro_commune;
	float biomasse_pechee_herb_commune_annuel;
	float biomasse_pechee_carniv_commune_annuel;
	float biomasse_pechee_herb_commune_annuel_r;
	float biomasse_pechee_carniv_commune_annuel_r;
	int nb_conflits_nuit_commune_annuel;
	int nb_conflits_jour_pecheur_commune_annuel;
	int nb_conflits_jour_tourisme_commune_annuel;
	int nb_quota_atteint_commune_annuel;
	
}


experiment exploration type: batch until: temps = temps_simulation repeat: 4 keep_seed: true{
	parameter debut_perturbation var: debut_perturbation among:[1000.0, 100000.0];
	parameter quota var:quota among: [0.5, 10000.0];
	parameter aide_financiere var:aide_financiere among: [0.0, 0.3];
	parameter mode_batch var:mode_batch <- true;
}


experiment complete_world type: gui {
	
	init {
		//create simulation with: [perturbation :: true, output_carte_lagon :: "../Outputs/Statu quo/COTS/Lagon/statu_quo_perturbation_lagon_", output_carte_commune :: "../Outputs/Statu quo/COTS/Communes/statu_quo_perturbation_commune_", output_csv_ecologique :: "../Outputs/Statu quo/COTS/statu_quo_perturbation_ecologique.csv", output_csv_cell :: "../Outputs/Statu quo/COTS/Cellules/statu_quo_perturbation_cell_", output_csv_pecheur :: "../Outputs/Statu quo/COTS/Pecheurs/statu_quo_perturbation_pecheur_"];
        
    }
	
	output {
				
		display "Chart corails - algues" refresh: every(14.0) {
			chart "Lagon" type: series size: {1.0, 0.5} {
				data "Corail" value: C_moyenne_lagon; 
				data "Turf" value: T_moyenne_lagon;				
			}			
			chart "Pente externe" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
				data "Corail" value: C_moyenne_pente_externe; 
				data "Turf" value: T_moyenne_pente_externe;
			}
		}
		
		display "Chart Pêche" refresh: every(1.0) {						
			chart "Biomasse de poissons pêchée" type: series {				
					data "Biomasse herbivores pente externe"  value: biomasse_pechee_herb_pente_externe color: #green marker_shape: marker_circle;
					data "Biomasse herbivores lagon"  value: biomasse_pechee_herb_lagon color: #green marker_shape: marker_square;					
					data "Biomasse carnivores pente externe"  value: biomasse_pechee_carniv_pente_externe color: #red marker_shape: marker_circle;
					data "Biomasse carnivores lagon"  value: biomasse_pechee_carniv_lagon color: #red marker_shape: marker_square;		
			}	
		}
		
		display "Chart Pêcheurs" refresh: every(1.0) {
			chart "Quotas atteint" type: series size: {1.0, 0.5}  {				
					data "Quotas atteint"  value: quota_atteint color: #black;
			}
			chart "Nb conflits lagon" type: series size: {1.0, 0.5} position: {0.0, 0.5} {				
					data "Nombre de conflits pêcheur de jour" value: conflits_jour_pecheur color: #orange;
					data "Nombre de conflits tourisme de jour" value: conflits_jour_tourisme color: #orange; 
					data "Nombre de conflits la nuit" value: conflits_nuit color: #black;
 			}	
		}
		
		display "Chart fishes" refresh: every(1.0) {
			chart "Biomasse de poissons lagon" type: series size: {1.0, 0.5} {				
					data "Herbivores" value: H_moyenne_lagon color: #green; 
					data "Corallivores" value: CO_moyenne_lagon color: #blue; 
					data "Carnivores" value: P_moyenne_lagon color: #red;					
			}			
			chart "Biomasse de poissons pente externe" type: series size: {1.0, 0.5} position: {0.0, 0.5} {				
					data "Herbivores" value: H_moyenne_pente_externe color: #green; 
					data "Corallivores" value: CO_moyenne_pente_externe color: #blue; 
					data "Carnivores" value: P_moyenne_pente_externe color: #red;		
			}
		}
		display "Chart fishes AMP" refresh: every(1.0) {
			chart "Biomasse de poissons AMP" type: series size: {1.0, 0.5} {				
					data "Herbivores" value: H_moyenne_AMP color: #green; 
					data "Corallivores" value: CO_moyenne_AMP color: #blue; 
					data "Carnivores" value: P_moyenne_AMP color: #red;					
			}			
			chart "Biomasse de poissons hors AMP" type: series size: {1.0, 0.5} position: {0.0, 0.5} {				
					data "Herbivores" value: H_moyenne_not_AMP color: #green; 
					data "Corallivores" value: CO_moyenne_not_AMP color: #blue; 
					data "Carnivores" value: P_moyenne_not_AMP color: #red;		
			}
		}
	}
}

experiment clipped_world type: gui {
	parameter clipped_data var:clipped_data <- true;
	
	init {
		//create simulation with:[biomasse_pechee_moyenne :: 0.8, output_carte_lagon :: "../Outputs/Statut quo/COTS/statut_quo_perturbation_lagon.shp", output_csv_district :: "../Outputs/Statut quo/COTS/statut_quo_perturbation_district.csv", output_txt :: "../Outputs/Statut quo/COTS/statut_quo_perturbation.txt", output_csv_ecologique :: "../Outputs/Statut quo/COTS/statut_quo_perturbation_ecologique.csv", output_csv_cell :: "../Outputs/Statut quo/COTS/statut_quo_perturbation_cell.csv", output_csv_pecheur :: "../Outputs/Statut quo/COTS/statut_quo_perturbation_pecheur.csv"];
        
    }
	
	output {
		
		display "Chart fishes" refresh: every(1.0) {
			chart "Biomasse de poissons lagon" type: series size: {1.0, 0.5} {				
					data "Herbivores" value: H_moyenne_lagon color: #green; 
					data "Corallivores" value: CO_moyenne_lagon color: #blue; 
					data "Carnivores" value: P_moyenne_lagon color: #red;					
			}			
			chart "Biomasse de poissons pente externe" type: series size: {1.0, 0.5} position: {0.0, 0.5} {				
					data "Herbivores" value: H_moyenne_pente_externe color: #green; 
					data "Corallivores" value: CO_moyenne_pente_externe color: #blue; 
					data "Carnivores" value: P_moyenne_pente_externe color: #red;		
			}
		}
		display "Chart fishes AMP" refresh: every(1.0) {
			chart "Biomasse de poissons AMP" type: series size: {1.0, 0.5} {				
					data "Herbivores" value: H_moyenne_AMP color: #green; 
					data "Corallivores" value: CO_moyenne_AMP color: #blue; 
					data "Carnivores" value: P_moyenne_AMP color: #red;					
			}			
			chart "Biomasse de poissons hors AMP" type: series size: {1.0, 0.5} position: {0.0, 0.5} {				
					data "Herbivores" value: H_moyenne_not_AMP color: #green; 
					data "Corallivores" value: CO_moyenne_not_AMP color: #blue; 
					data "Carnivores" value: P_moyenne_not_AMP color: #red;		
			}
		}
		
		
		display "Chart corails - algues" refresh: every(1.0) {
			chart "Lagon" type: series size: {1.0, 0.5} {
				data "Corail" value: C_moyenne_lagon; 
				data "Turf" value: T_moyenne_lagon;				
			}			
			chart "Pente externe" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
				data "Corail" value: C_moyenne_pente_externe; 
				data "Turf" value: T_moyenne_pente_externe;
			}
		}
		
		display "Chart Pêche" refresh: every(1.0) {						
			chart "Biomasse de poissons pêchée" type: series {				
					data "Biomasse herbivores pente externe"  value: biomasse_pechee_herb_pente_externe color: #green marker_shape: marker_circle;
					data "Biomasse herbivores lagon"  value: biomasse_pechee_herb_lagon color: #green marker_shape: marker_square;					
					data "Biomasse carnivores pente externe"  value: biomasse_pechee_carniv_pente_externe color: #red marker_shape: marker_circle;
					data "Biomasse carnivores lagon"  value: biomasse_pechee_carniv_lagon color: #red marker_shape: marker_square;		
			}	
		}
		
		display "Chart Pêcheurs" refresh: every(1.0) {
			chart "Quotas atteint" type: series size: {1.0, 0.5}  {				
					data "Quotas atteint"  value: quota_atteint color: #black;
			}
			chart "Nb conflits lagon" type: series size: {1.0, 0.5} position: {0.0, 0.5} {
					data "Nombre de conflits la nuit lagon" value: conflits_nuit color: #black;			
					data "Nombre de conflits pêcheur de jour lagon" value: conflits_jour_pecheur color: #orange;
					data "Nombre de conflits tourisme de jour lagon" value: conflits_jour_tourisme color: #orange; 
 			}	
		}
		
	}
}