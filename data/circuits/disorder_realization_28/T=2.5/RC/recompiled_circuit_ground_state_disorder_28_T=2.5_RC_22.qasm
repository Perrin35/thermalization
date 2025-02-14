OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82792038) q[0];
sx q[0];
rz(-1.6452687) q[0];
sx q[0];
rz(-0.21284719) q[0];
rz(-1.4053474) q[1];
sx q[1];
rz(-1.5264629) q[1];
sx q[1];
rz(-0.57797536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5151857) q[0];
sx q[0];
rz(-0.66045633) q[0];
sx q[0];
rz(-0.68645262) q[0];
rz(1.3499267) q[2];
sx q[2];
rz(-0.71227461) q[2];
sx q[2];
rz(0.28755915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2541952) q[1];
sx q[1];
rz(-2.7294558) q[1];
sx q[1];
rz(-1.5855476) q[1];
x q[2];
rz(-2.7733581) q[3];
sx q[3];
rz(-1.7126178) q[3];
sx q[3];
rz(-0.54556812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9389682) q[2];
sx q[2];
rz(-1.5105931) q[2];
sx q[2];
rz(0.6187588) q[2];
rz(0.42801157) q[3];
sx q[3];
rz(-1.2330387) q[3];
sx q[3];
rz(1.7091883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0851058) q[0];
sx q[0];
rz(-3.1054057) q[0];
sx q[0];
rz(-2.4888743) q[0];
rz(2.808049) q[1];
sx q[1];
rz(-0.8077375) q[1];
sx q[1];
rz(2.6127167) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655114) q[0];
sx q[0];
rz(-1.2032857) q[0];
sx q[0];
rz(-1.0239081) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0300601) q[2];
sx q[2];
rz(-2.1443488) q[2];
sx q[2];
rz(-2.3071764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6436504) q[1];
sx q[1];
rz(-1.4406246) q[1];
sx q[1];
rz(-1.4243717) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5020735) q[3];
sx q[3];
rz(-2.4009224) q[3];
sx q[3];
rz(3.1225151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35398206) q[2];
sx q[2];
rz(-1.3926316) q[2];
sx q[2];
rz(-0.37484136) q[2];
rz(1.2104872) q[3];
sx q[3];
rz(-2.6755302) q[3];
sx q[3];
rz(-1.2208285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0354563) q[0];
sx q[0];
rz(-1.184967) q[0];
sx q[0];
rz(2.5842066) q[0];
rz(-0.71929559) q[1];
sx q[1];
rz(-2.6917515) q[1];
sx q[1];
rz(-0.61190277) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9111262) q[0];
sx q[0];
rz(-1.9042) q[0];
sx q[0];
rz(0.35349288) q[0];
x q[1];
rz(-0.8863897) q[2];
sx q[2];
rz(-1.4006613) q[2];
sx q[2];
rz(-1.1332133) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53781834) q[1];
sx q[1];
rz(-2.1421931) q[1];
sx q[1];
rz(-0.091544108) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40308904) q[3];
sx q[3];
rz(-1.7397907) q[3];
sx q[3];
rz(2.4141261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3376075) q[2];
sx q[2];
rz(-1.7943725) q[2];
sx q[2];
rz(0.95532974) q[2];
rz(2.5721926) q[3];
sx q[3];
rz(-1.1207213) q[3];
sx q[3];
rz(-1.3620522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78205458) q[0];
sx q[0];
rz(-2.8466917) q[0];
sx q[0];
rz(-3.0468347) q[0];
rz(1.5358198) q[1];
sx q[1];
rz(-1.9784119) q[1];
sx q[1];
rz(0.43704978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9511918) q[0];
sx q[0];
rz(-0.13057183) q[0];
sx q[0];
rz(-0.87223069) q[0];
rz(-2.7669698) q[2];
sx q[2];
rz(-0.46430507) q[2];
sx q[2];
rz(1.855206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0298095) q[1];
sx q[1];
rz(-2.8800721) q[1];
sx q[1];
rz(-1.5716121) q[1];
x q[2];
rz(-2.0548204) q[3];
sx q[3];
rz(-0.71992249) q[3];
sx q[3];
rz(-2.7663224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.36865083) q[2];
sx q[2];
rz(-2.3029885) q[2];
sx q[2];
rz(1.3402026) q[2];
rz(-2.7784427) q[3];
sx q[3];
rz(-2.4996417) q[3];
sx q[3];
rz(-0.44410458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7813613) q[0];
sx q[0];
rz(-2.4960127) q[0];
sx q[0];
rz(3.0907104) q[0];
rz(-2.6119192) q[1];
sx q[1];
rz(-2.5310204) q[1];
sx q[1];
rz(3.1231336) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7537024) q[0];
sx q[0];
rz(-1.876693) q[0];
sx q[0];
rz(1.0359156) q[0];
rz(-pi) q[1];
rz(0.71645488) q[2];
sx q[2];
rz(-1.418651) q[2];
sx q[2];
rz(-2.1670451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27996846) q[1];
sx q[1];
rz(-1.8496152) q[1];
sx q[1];
rz(-2.2303225) q[1];
x q[2];
rz(-2.6191447) q[3];
sx q[3];
rz(-1.4178489) q[3];
sx q[3];
rz(0.82081036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.20794836) q[2];
sx q[2];
rz(-1.7174145) q[2];
sx q[2];
rz(-1.2733744) q[2];
rz(-0.4452855) q[3];
sx q[3];
rz(-0.71355009) q[3];
sx q[3];
rz(2.2860897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2263366) q[0];
sx q[0];
rz(-2.7610918) q[0];
sx q[0];
rz(-0.29391995) q[0];
rz(1.2433012) q[1];
sx q[1];
rz(-1.2970122) q[1];
sx q[1];
rz(1.4909202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7939756) q[0];
sx q[0];
rz(-1.4039984) q[0];
sx q[0];
rz(2.1991437) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8155509) q[2];
sx q[2];
rz(-1.4801822) q[2];
sx q[2];
rz(-0.88256529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6676644) q[1];
sx q[1];
rz(-1.7662303) q[1];
sx q[1];
rz(0.21904314) q[1];
rz(-pi) q[2];
rz(-1.1209247) q[3];
sx q[3];
rz(-0.48512019) q[3];
sx q[3];
rz(-2.7435722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9320429) q[2];
sx q[2];
rz(-1.2806634) q[2];
sx q[2];
rz(-3.0397084) q[2];
rz(-2.3515676) q[3];
sx q[3];
rz(-0.70370379) q[3];
sx q[3];
rz(-1.8858645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0027960221) q[0];
sx q[0];
rz(-0.094745435) q[0];
sx q[0];
rz(2.5139659) q[0];
rz(-1.3603285) q[1];
sx q[1];
rz(-1.5831213) q[1];
sx q[1];
rz(2.9291709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14721641) q[0];
sx q[0];
rz(-1.205814) q[0];
sx q[0];
rz(2.2969691) q[0];
x q[1];
rz(-0.27842267) q[2];
sx q[2];
rz(-2.1489193) q[2];
sx q[2];
rz(-2.9591564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46225658) q[1];
sx q[1];
rz(-1.0733593) q[1];
sx q[1];
rz(-0.21381198) q[1];
x q[2];
rz(-1.095781) q[3];
sx q[3];
rz(-1.1355564) q[3];
sx q[3];
rz(2.4839885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75954306) q[2];
sx q[2];
rz(-2.2727649) q[2];
sx q[2];
rz(-2.7810435) q[2];
rz(-0.77066317) q[3];
sx q[3];
rz(-1.9375075) q[3];
sx q[3];
rz(-2.8785021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78779546) q[0];
sx q[0];
rz(-2.0833092) q[0];
sx q[0];
rz(-1.298792) q[0];
rz(-0.53954387) q[1];
sx q[1];
rz(-1.455247) q[1];
sx q[1];
rz(-1.4849327) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4277219) q[0];
sx q[0];
rz(-1.026928) q[0];
sx q[0];
rz(-1.1975678) q[0];
x q[1];
rz(-2.7168112) q[2];
sx q[2];
rz(-1.3662587) q[2];
sx q[2];
rz(1.0077602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.114552) q[1];
sx q[1];
rz(-1.5590347) q[1];
sx q[1];
rz(-1.8742754) q[1];
rz(2.4267107) q[3];
sx q[3];
rz(-1.2586888) q[3];
sx q[3];
rz(-1.1955137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7654884) q[2];
sx q[2];
rz(-1.8643943) q[2];
sx q[2];
rz(-2.6831324) q[2];
rz(1.9930528) q[3];
sx q[3];
rz(-0.12980041) q[3];
sx q[3];
rz(-0.54500088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.604436) q[0];
sx q[0];
rz(-0.24188365) q[0];
sx q[0];
rz(2.3271374) q[0];
rz(2.3711329) q[1];
sx q[1];
rz(-1.6043681) q[1];
sx q[1];
rz(1.3260215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6961937) q[0];
sx q[0];
rz(-2.8210381) q[0];
sx q[0];
rz(-2.205109) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4374323) q[2];
sx q[2];
rz(-1.6495145) q[2];
sx q[2];
rz(0.12264473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5863813) q[1];
sx q[1];
rz(-1.5813236) q[1];
sx q[1];
rz(1.7278683) q[1];
x q[2];
rz(1.7538449) q[3];
sx q[3];
rz(-1.7244851) q[3];
sx q[3];
rz(-2.1240687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.068437964) q[2];
sx q[2];
rz(-2.1775553) q[2];
sx q[2];
rz(-0.55802074) q[2];
rz(2.8716715) q[3];
sx q[3];
rz(-0.25756613) q[3];
sx q[3];
rz(-2.423438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.54616) q[0];
sx q[0];
rz(-1.8159001) q[0];
sx q[0];
rz(-0.85920715) q[0];
rz(-2.4190306) q[1];
sx q[1];
rz(-0.95046202) q[1];
sx q[1];
rz(1.7019466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50558907) q[0];
sx q[0];
rz(-3.1312864) q[0];
sx q[0];
rz(2.206102) q[0];
x q[1];
rz(1.7990803) q[2];
sx q[2];
rz(-1.205027) q[2];
sx q[2];
rz(-1.7342784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0175914) q[1];
sx q[1];
rz(-0.81087001) q[1];
sx q[1];
rz(-1.3254741) q[1];
rz(-0.4018114) q[3];
sx q[3];
rz(-1.7932745) q[3];
sx q[3];
rz(-2.1232186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47448164) q[2];
sx q[2];
rz(-1.8626532) q[2];
sx q[2];
rz(-1.4943592) q[2];
rz(-0.4161559) q[3];
sx q[3];
rz(-1.6437203) q[3];
sx q[3];
rz(-0.070629899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304929) q[0];
sx q[0];
rz(-2.6051705) q[0];
sx q[0];
rz(2.312533) q[0];
rz(-0.88498712) q[1];
sx q[1];
rz(-0.61274715) q[1];
sx q[1];
rz(-1.3428584) q[1];
rz(-1.5985684) q[2];
sx q[2];
rz(-2.118961) q[2];
sx q[2];
rz(-1.9365666) q[2];
rz(2.3571499) q[3];
sx q[3];
rz(-1.389849) q[3];
sx q[3];
rz(0.72985284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
