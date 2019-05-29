--------------------------------------------------------
--
-- Module: General
--
--------------------------------------------------------
	
REVOKE ALL ON tableinfo FROM PUBLIC;
REVOKE ALL ON db FROM PUBLIC;
REVOKE ALL ON dbxref FROM PUBLIC;
REVOKE ALL ON project FROM PUBLIC;

--------------------------------------------------------
--
-- Module: CV
--
--------------------------------------------------------
REVOKE ALL ON cv FROM PUBLIC;
REVOKE ALL ON cvterm FROM PUBLIC;
REVOKE ALL ON cvterm_relationship FROM PUBLIC;
REVOKE ALL ON cvtermpath FROM PUBLIC;
REVOKE ALL ON cvtermsynonym FROM PUBLIC;
REVOKE ALL ON cvterm_dbxref FROM PUBLIC;
REVOKE ALL ON cvtermprop FROM PUBLIC;
REVOKE ALL ON dbxrefprop FROM PUBLIC;

--------------------------------------------------------
--
-- Module: Pub
--
--------------------------------------------------------
REVOKE ALL ON pub FROM PUBLIC;
REVOKE ALL ON pub_relationship FROM PUBLIC;
REVOKE ALL ON pub_dbxref FROM PUBLIC;
REVOKE ALL ON pubauthor FROM PUBLIC;
REVOKE ALL ON pubprop FROM PUBLIC;

--------------------------------------------------------
--
-- Module: Organism
--
--------------------------------------------------------
REVOKE ALL ON organism FROM PUBLIC;
REVOKE ALL ON organism_dbxref FROM PUBLIC;
REVOKE ALL ON organismprop FROM PUBLIC;

--------------------------------------------------------
--
-- Module: Sequence
--
--------------------------------------------------------
REVOKE ALL ON feature FROM PUBLIC;
REVOKE ALL ON featureloc FROM PUBLIC;
REVOKE ALL ON feature_pub FROM PUBLIC;
REVOKE ALL ON featureprop FROM PUBLIC;
REVOKE ALL ON featureprop_pub FROM PUBLIC;
REVOKE ALL ON feature_dbxref FROM PUBLIC;
REVOKE ALL ON feature_relationship FROM PUBLIC;
REVOKE ALL ON feature_relationship_pub FROM PUBLIC;
REVOKE ALL ON feature_relationshipprop FROM PUBLIC;
REVOKE ALL ON feature_relprop_pub FROM PUBLIC;
REVOKE ALL ON feature_cvterm FROM PUBLIC;
REVOKE ALL ON feature_cvtermprop FROM PUBLIC;
REVOKE ALL ON feature_cvterm_dbxref FROM PUBLIC;
REVOKE ALL ON feature_cvterm_pub FROM PUBLIC;
REVOKE ALL ON synonym FROM PUBLIC;
REVOKE ALL ON feature_synonym FROM PUBLIC;

--------------------------------------------------------
--
-- Module: Computational Analysis
--
--------------------------------------------------------
REVOKE ALL ON analysis FROM PUBLIC;
REVOKE ALL ON analysisprop FROM PUBLIC;
REVOKE ALL ON analysisfeature FROM PUBLIC;



--------------------------------------------------------
--
-- Module: Phylogeny
--
--------------------------------------------------------

REVOKE ALL ON phylotree FROM PUBLIC;
REVOKE ALL ON phylotree_pub FROM PUBLIC;
REVOKE ALL ON phylonode FROM PUBLIC;
REVOKE ALL ON phylonode_dbxref FROM PUBLIC;
REVOKE ALL ON phylonode_pub FROM PUBLIC;
REVOKE ALL ON phylonode_organism FROM PUBLIC;
REVOKE ALL ON phylonodeprop FROM PUBLIC;
REVOKE ALL ON phylonode_relationship FROM PUBLIC;

-------------------------------------------------------
--
-- Module: Infection
--
-------------------------------------------------------
REVOKE ALL ON infection FROM PUBLIC;
REVOKE ALL ON infectionprop FROM PUBLIC;
REVOKE ALL ON infection_cvterm FROM PUBLIC;
REVOKE ALL ON infection_dbxref FROM PUBLIC;
REVOKE ALL ON transmission FROM PUBLIC;
REVOKE ALL ON transmissionprop FROM PUBLIC;
REVOKE ALL ON transmission_cvterm FROM PUBLIC;
REVOKE ALL ON transmission_dbxref FROM PUBLIC;
REVOKE ALL ON incident FROM PUBLIC;
REVOKE ALL ON incidentprop FROM PUBLIC;
REVOKE ALL ON incident_cvterm FROM PUBLIC;
REVOKE ALL ON incident_dbxref FROM PUBLIC;
REVOKE ALL ON incident_relationship FROM PUBLIC;
REVOKE ALL ON infection_pub FROM PUBLIC;
REVOKE ALL ON transmission_pub FROM PUBLIC;
REVOKE ALL ON incident_pub FROM PUBLIC;

------------------------------------------------------
--
-- Module: Contact
--
------------------------------------------------------
REVOKE ALL ON contact FROM PUBLIC;
REVOKE ALL ON contactprop FROM PUBLIC;
REVOKE ALL ON contact_relationship FROM PUBLIC;


------------------------------------------------------
--
-- Module: Genotype
--
------------------------------------------------------
REVOKE ALL ON genotype FROM PUBLIC;
REVOKE ALL ON feature_genotype FROM PUBLIC;
REVOKE ALL ON environment FROM PUBLIC;
REVOKE ALL ON environment_cvterm FROM PUBLIC;
REVOKE ALL ON phenstatement FROM PUBLIC;
REVOKE ALL ON phendesc FROM PUBLIC;
REVOKE ALL ON phenotype_comparison FROM PUBLIC;
REVOKE ALL ON phenotype_comparison_cvterm FROM PUBLIC;


-----------------------------------------------------
--
-- Module: Stock
--
-----------------------------------------------------
REVOKE ALL ON stock FROM PUBLIC;
REVOKE ALL ON stock_pub FROM PUBLIC;
REVOKE ALL ON stockprop FROM PUBLIC;
REVOKE ALL ON stockprop_pub FROM PUBLIC;
REVOKE ALL ON stock_relationship FROM PUBLIC;
REVOKE ALL ON stock_relationship_pub FROM PUBLIC;
REVOKE ALL ON stock_dbxref FROM PUBLIC;
REVOKE ALL ON stock_cvterm FROM PUBLIC;
REVOKE ALL ON stock_genotype FROM PUBLIC;
REVOKE ALL ON stockcollection FROM PUBLIC;
REVOKE ALL ON stockcollectionprop FROM PUBLIC;

