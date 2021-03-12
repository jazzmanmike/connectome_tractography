# connectome_tractography
 Matlab code for connectome analysis with tractography data
 
 See the help header in Matlab for full details. 
 
 Dependencies include: BCT 2019_03_03, contest, versatility, matlab_bgl, powerlaws_full, schemaball, BrainNetViewer
 
 Inputs: connectivity matrix & parcellation template co-ordinates
 
 To get the above either use tract_van.sh for processing of tractrography data or alternatively just use your own data. 
 
 To start just need to enter the patient path in section A1.

            %% A1. Load data

            %This should be the only part required to be set manually

            %Directory
            directory = '/path_to_my_data';

            %Patient ID
            patientID = 'my_patient_ID';

            %Template name
            template = 'AAL90';
