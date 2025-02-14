OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0877617) q[0];
sx q[0];
rz(-0.55019903) q[0];
sx q[0];
rz(1.7734111) q[0];
rz(2.6693681) q[1];
sx q[1];
rz(-0.11657403) q[1];
sx q[1];
rz(0.20224686) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7769788) q[0];
sx q[0];
rz(-1.3585416) q[0];
sx q[0];
rz(-1.3670761) q[0];
rz(-pi) q[1];
rz(-0.80533452) q[2];
sx q[2];
rz(-0.67395681) q[2];
sx q[2];
rz(1.0588624) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.80965197) q[1];
sx q[1];
rz(-0.9729079) q[1];
sx q[1];
rz(0.34418587) q[1];
rz(-pi) q[2];
rz(-1.4957761) q[3];
sx q[3];
rz(-2.6743691) q[3];
sx q[3];
rz(-2.644616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54174417) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8180654) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(0.7793119) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(1.13387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72778217) q[0];
sx q[0];
rz(-1.4725794) q[0];
sx q[0];
rz(-2.006524) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0236565) q[2];
sx q[2];
rz(-1.6291766) q[2];
sx q[2];
rz(-2.749325) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6821377) q[1];
sx q[1];
rz(-1.3054779) q[1];
sx q[1];
rz(3.0299761) q[1];
x q[2];
rz(2.0800703) q[3];
sx q[3];
rz(-1.5699374) q[3];
sx q[3];
rz(-2.1775376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69089943) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(-2.9449985) q[2];
rz(2.172566) q[3];
sx q[3];
rz(-1.5195547) q[3];
sx q[3];
rz(-1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3900688) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(2.4098136) q[0];
rz(-0.36830184) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-0.36645737) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5598504) q[0];
sx q[0];
rz(-1.7479715) q[0];
sx q[0];
rz(0.26818256) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7973034) q[2];
sx q[2];
rz(-1.3346121) q[2];
sx q[2];
rz(0.073848595) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11135396) q[1];
sx q[1];
rz(-2.0890279) q[1];
sx q[1];
rz(0.49211647) q[1];
rz(-0.89160925) q[3];
sx q[3];
rz(-2.3162033) q[3];
sx q[3];
rz(2.2872383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9590108) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(0.89243531) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733805) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(0.40661231) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(2.3060395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0428187) q[0];
sx q[0];
rz(-1.6102428) q[0];
sx q[0];
rz(-1.7516319) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0876483) q[2];
sx q[2];
rz(-2.7479541) q[2];
sx q[2];
rz(-2.7700499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16627993) q[1];
sx q[1];
rz(-2.0034092) q[1];
sx q[1];
rz(1.4090502) q[1];
rz(-pi) q[2];
rz(1.1855679) q[3];
sx q[3];
rz(-1.480417) q[3];
sx q[3];
rz(2.4848695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(0.92140222) q[2];
rz(-1.7737927) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88857404) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(2.9827523) q[0];
rz(1.6244434) q[1];
sx q[1];
rz(-2.8327063) q[1];
sx q[1];
rz(-1.812017) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3811262) q[0];
sx q[0];
rz(-1.5601349) q[0];
sx q[0];
rz(0.018254781) q[0];
rz(-pi) q[1];
rz(0.23933584) q[2];
sx q[2];
rz(-2.1194601) q[2];
sx q[2];
rz(1.2557097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.98011866) q[1];
sx q[1];
rz(-1.1921645) q[1];
sx q[1];
rz(2.9254854) q[1];
x q[2];
rz(1.8908126) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(-1.0726014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57881957) q[2];
sx q[2];
rz(-2.4476738) q[2];
sx q[2];
rz(-2.5895183) q[2];
rz(2.3479346) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8296705) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(2.435834) q[0];
rz(0.13748473) q[1];
sx q[1];
rz(-2.4973713) q[1];
sx q[1];
rz(-2.4286043) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.997809) q[0];
sx q[0];
rz(-1.1857872) q[0];
sx q[0];
rz(1.817817) q[0];
rz(2.627264) q[2];
sx q[2];
rz(-1.1900969) q[2];
sx q[2];
rz(-1.9388388) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29999816) q[1];
sx q[1];
rz(-1.9683241) q[1];
sx q[1];
rz(0.10511153) q[1];
x q[2];
rz(1.2014404) q[3];
sx q[3];
rz(-0.31285646) q[3];
sx q[3];
rz(0.85858708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2116427) q[2];
sx q[2];
rz(-2.2723618) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(-0.12187135) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(-1.6596644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26509869) q[0];
sx q[0];
rz(-2.5229186) q[0];
sx q[0];
rz(-2.4342243) q[0];
rz(-2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.3622989) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4934568) q[0];
sx q[0];
rz(-1.8168983) q[0];
sx q[0];
rz(1.2854693) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84378924) q[2];
sx q[2];
rz(-1.8966881) q[2];
sx q[2];
rz(-1.0032723) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41304663) q[1];
sx q[1];
rz(-0.73592454) q[1];
sx q[1];
rz(2.3152183) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45622042) q[3];
sx q[3];
rz(-1.4672605) q[3];
sx q[3];
rz(0.52952784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97544396) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(0.56900209) q[2];
rz(2.1042018) q[3];
sx q[3];
rz(-2.2834957) q[3];
sx q[3];
rz(-2.67498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55050945) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(-1.9367223) q[0];
rz(2.6288746) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(-2.9606294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2406684) q[0];
sx q[0];
rz(-2.5818733) q[0];
sx q[0];
rz(0.75738917) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7920807) q[2];
sx q[2];
rz(-1.4415359) q[2];
sx q[2];
rz(-1.845128) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6540203) q[1];
sx q[1];
rz(-2.1250238) q[1];
sx q[1];
rz(0.2481064) q[1];
rz(-pi) q[2];
rz(0.11060235) q[3];
sx q[3];
rz(-1.6440653) q[3];
sx q[3];
rz(-1.4721118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60244954) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(-3.0446206) q[2];
rz(2.5868331) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(-0.35840148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612514) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(0.10777792) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(-1.617618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.658889) q[0];
sx q[0];
rz(-1.3388757) q[0];
sx q[0];
rz(-0.13302444) q[0];
x q[1];
rz(-0.6438821) q[2];
sx q[2];
rz(-1.9877627) q[2];
sx q[2];
rz(1.3326926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7858582) q[1];
sx q[1];
rz(-1.5275035) q[1];
sx q[1];
rz(1.1086575) q[1];
rz(1.518256) q[3];
sx q[3];
rz(-1.0797622) q[3];
sx q[3];
rz(2.4201916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9109351) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-2.0111734) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(-2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0216574) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(0.94104952) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(-1.2233268) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4775765) q[0];
sx q[0];
rz(-1.7682791) q[0];
sx q[0];
rz(-0.5911747) q[0];
rz(-2.1603196) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(-1.5379932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0666484) q[1];
sx q[1];
rz(-2.7171591) q[1];
sx q[1];
rz(-0.17129274) q[1];
rz(-1.2480601) q[3];
sx q[3];
rz(-0.56005037) q[3];
sx q[3];
rz(-1.3572451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90784043) q[2];
sx q[2];
rz(-1.7641822) q[2];
sx q[2];
rz(-0.09093786) q[2];
rz(0.20445538) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37888708) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(1.0574404) q[2];
sx q[2];
rz(-0.79887894) q[2];
sx q[2];
rz(0.32499921) q[2];
rz(2.8818535) q[3];
sx q[3];
rz(-2.9208214) q[3];
sx q[3];
rz(-0.95820273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
