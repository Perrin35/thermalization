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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(0.51796335) q[1];
sx q[1];
rz(-1.0022751) q[1];
sx q[1];
rz(0.60751539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54017055) q[0];
sx q[0];
rz(-0.54781944) q[0];
sx q[0];
rz(2.026985) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0614228) q[2];
sx q[2];
rz(-1.4685838) q[2];
sx q[2];
rz(-1.4776023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3594397) q[1];
sx q[1];
rz(-0.50737587) q[1];
sx q[1];
rz(-0.33276593) q[1];
x q[2];
rz(2.257454) q[3];
sx q[3];
rz(-2.5885785) q[3];
sx q[3];
rz(0.64698863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6174378) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(0.080862008) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2422159) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(-2.2826165) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(-1.7346409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1340356) q[0];
sx q[0];
rz(-1.3430183) q[0];
sx q[0];
rz(-2.0640949) q[0];
x q[1];
rz(1.2386049) q[2];
sx q[2];
rz(-1.6055577) q[2];
sx q[2];
rz(-1.3369651) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9216361) q[1];
sx q[1];
rz(-1.0590648) q[1];
sx q[1];
rz(3.0380766) q[1];
rz(2.5646016) q[3];
sx q[3];
rz(-1.2170047) q[3];
sx q[3];
rz(2.4959223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0824288) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.814369) q[2];
rz(-2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(-0.41675848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5752207) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(2.7161993) q[0];
rz(1.7644024) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4017556) q[0];
sx q[0];
rz(-2.4889737) q[0];
sx q[0];
rz(-2.9463861) q[0];
rz(0.70830958) q[2];
sx q[2];
rz(-1.5454486) q[2];
sx q[2];
rz(0.011653221) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9454398) q[1];
sx q[1];
rz(-1.1588133) q[1];
sx q[1];
rz(-0.71228551) q[1];
x q[2];
rz(1.9886964) q[3];
sx q[3];
rz(-0.66158453) q[3];
sx q[3];
rz(2.083287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(1.8505992) q[2];
rz(1.7839606) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(-1.4972081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62035471) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(-2.3714016) q[0];
rz(2.1417446) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(1.3410478) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4907341) q[0];
sx q[0];
rz(-0.58877173) q[0];
sx q[0];
rz(0.38932271) q[0];
rz(-pi) q[1];
rz(0.88444986) q[2];
sx q[2];
rz(-1.2053066) q[2];
sx q[2];
rz(-2.3523503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.136678) q[1];
sx q[1];
rz(-1.305849) q[1];
sx q[1];
rz(-2.8317004) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2433231) q[3];
sx q[3];
rz(-0.30042111) q[3];
sx q[3];
rz(1.6302366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(2.3731903) q[2];
rz(3.0411804) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860344) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.7598033) q[1];
sx q[1];
rz(-0.58265668) q[1];
sx q[1];
rz(1.7452128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5982626) q[0];
sx q[0];
rz(-0.51307438) q[0];
sx q[0];
rz(-2.3725879) q[0];
x q[1];
rz(-2.7740741) q[2];
sx q[2];
rz(-0.43635363) q[2];
sx q[2];
rz(0.23772109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5310881) q[1];
sx q[1];
rz(-0.63734326) q[1];
sx q[1];
rz(0.12048851) q[1];
rz(-pi) q[2];
rz(3.1110686) q[3];
sx q[3];
rz(-0.82477335) q[3];
sx q[3];
rz(-1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22623006) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(1.3876312) q[3];
sx q[3];
rz(-0.97532719) q[3];
sx q[3];
rz(1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607894) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(1.3355108) q[0];
rz(-0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(2.9551771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53509287) q[0];
sx q[0];
rz(-2.4267174) q[0];
sx q[0];
rz(2.2660648) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30419402) q[2];
sx q[2];
rz(-0.59202164) q[2];
sx q[2];
rz(-1.5409868) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5552206) q[1];
sx q[1];
rz(-2.456291) q[1];
sx q[1];
rz(1.9577115) q[1];
rz(-0.78808727) q[3];
sx q[3];
rz(-2.4895146) q[3];
sx q[3];
rz(-1.8502082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4621801) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.818044) q[2];
rz(1.9994252) q[3];
sx q[3];
rz(-2.3565632) q[3];
sx q[3];
rz(-0.89046684) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-2.9587726) q[0];
sx q[0];
rz(-0.63419813) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(2.1814836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8909559) q[0];
sx q[0];
rz(-1.5433558) q[0];
sx q[0];
rz(-1.9695884) q[0];
rz(0.30107408) q[2];
sx q[2];
rz(-1.5631526) q[2];
sx q[2];
rz(2.6507225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3030678) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(1.5477033) q[1];
rz(-pi) q[2];
rz(2.2957357) q[3];
sx q[3];
rz(-2.1921033) q[3];
sx q[3];
rz(1.814807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1977957) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(-1.5956399) q[2];
rz(-1.724285) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(-1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5328131) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(0.97484318) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625789) q[0];
sx q[0];
rz(-0.65555182) q[0];
sx q[0];
rz(-2.2961839) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5707827) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(0.9521614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(-1.8046677) q[1];
x q[2];
rz(-0.49174546) q[3];
sx q[3];
rz(-1.5881268) q[3];
sx q[3];
rz(2.8239856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10890266) q[2];
sx q[2];
rz(-0.87778512) q[2];
sx q[2];
rz(-0.6558134) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-0.18812215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.6498097) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(-2.6731491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830848) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(1.5751189) q[0];
rz(-0.55267398) q[2];
sx q[2];
rz(-0.50853339) q[2];
sx q[2];
rz(2.7486211) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2704308) q[1];
sx q[1];
rz(-0.79177815) q[1];
sx q[1];
rz(-1.3429848) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28161093) q[3];
sx q[3];
rz(-1.3002031) q[3];
sx q[3];
rz(0.56265807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.046772) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(2.2136733) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(2.0637276) q[0];
rz(0.36733356) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.2841388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3782702) q[0];
sx q[0];
rz(-2.5978932) q[0];
sx q[0];
rz(2.0732353) q[0];
x q[1];
rz(0.91725332) q[2];
sx q[2];
rz(-1.6628569) q[2];
sx q[2];
rz(1.9815418) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2312647) q[1];
sx q[1];
rz(-2.961425) q[1];
sx q[1];
rz(-0.50039165) q[1];
rz(-2.4068042) q[3];
sx q[3];
rz(-1.4126724) q[3];
sx q[3];
rz(1.2722335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1401356) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(-2.7122811) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(1.8163053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975288) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(2.5979832) q[2];
sx q[2];
rz(-1.3610441) q[2];
sx q[2];
rz(1.3112031) q[2];
rz(-2.024509) q[3];
sx q[3];
rz(-2.4632005) q[3];
sx q[3];
rz(1.7076422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
