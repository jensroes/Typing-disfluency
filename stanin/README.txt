# Brief description of all Stan models.

[1] "ARKMoGpptbgs.stan"                        
Autoregressive model with two mixture components with random effects for bigrams and participantsand autoregressor for each participant. All autoregression was first order. Mixing proportion was allowed to vary by participant.

[2] "ARKpptbgs.stan"                           
Autoregressive model with random effects for bigrams and participants and autoregressor for each participant. All autoregression was first order.

[3] "LMMbigramintercepts.stan"                 
Linear mixed-effects model with random effects for bigrams and participants.             

[4] "MoGpptsbigramintercepts.stan"            
Mixture model with two mixture components and random effects for bigrams and participants. 

[5] "MoGpptsbigraminterceptsdelta_betafix.stan"
Same as 4 but slowdown delta (second mixture component) was allowed to vary by participant instead of random intercepts.

[6] "MoGpptsbigraminterceptsdelta.stan"        
Same as 4 but slowdown delta (second mixture component) was allowed to vary by participant.

[7] "MoGpptsbigraminterceptsdelta2.stan" 
Same as 4 but slowdown delta (second mixture component) was allowed to vary by participant while mixing proportion theta was help constant.

