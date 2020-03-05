
# README

## OBS DATA
The file OBS_data.csv is a table with all the behavioral events registered by an observer, within a group of 20 Guinea baboons living in an enclosure of the Primate Center (CNRS-UPS846) in Rousset-sur-Arc (France), between June 13th 2019 and July 10th 2019.

The table consists of 8 columns:

- DateTime = Time stamp of the event, namely the moment the observed behavior was registered. In case of STATE events (events with duration > 0), it refers to the beginning of the behavior;

- Actor = The name of the actor;

- Recipient = The name of the individual the Actor is acting upon;

- Behavior = The behavior the Actor. 14 types of behaviors are registered:'Resting', 'Grooming', 'Presenting','Playing with', 'Grunting-Lipsmacking', 'Supplanting','Threatening', 'Submission', 'Touching', 'Avoiding', 'Attacking','Carrying', 'Embracing', 'Mounting', 'Copulating', 'Chasing'. In addition two other categories were included: 'Invisible' and 'Other';

- Category = The classification of the behavior. It can be 'Affiliative', 'Agonistic', 'Other'; 

- Duration = Duration of the observed behavior. POINT events have no duration;

- Localisation =  Zone of the enclosure where the observed behavior takes place; 

- Point = indicates if the event is a POINT event (YES) or a STATE event (NO).

## RFID DATA

The file RFID_data.csv is a table containing RFID contacts data recorded in the same period (June 13th 2019 - July 10th 2019).
The contacts are aggregated in time intervals of 20 seconds. They refer to couples of wearable proximity sensors (http://www.sociopatterns.org/), worn by 13 of the 20 individuals cited above. Data are registered by readers positioned around the enclosure. 
The name of the tags were replaced with the name of the individual wearing the tag, to be consistent with the observational data. 
The table consists of 4 columns: 

- t = time of the beginning of the contact in Epoch format (Unix timestamps);

- i = Name of the first individual;

- j = Name of the second individual; 

- DateTime = time converted in Time stamp format.
