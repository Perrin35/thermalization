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
rz(0.55819297) q[0];
sx q[0];
rz(3.8605122) q[0];
sx q[0];
rz(10.100848) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(2.2303384) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9118496) q[0];
sx q[0];
rz(-1.5930682) q[0];
sx q[0];
rz(-3.1215206) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5160271) q[2];
sx q[2];
rz(-1.5766597) q[2];
sx q[2];
rz(-1.0926343) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75685749) q[1];
sx q[1];
rz(-2.5674106) q[1];
sx q[1];
rz(3.0822445) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4449213) q[3];
sx q[3];
rz(-0.81842234) q[3];
sx q[3];
rz(0.5054121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91997826) q[2];
sx q[2];
rz(-0.45404926) q[2];
sx q[2];
rz(-0.81217074) q[2];
rz(3.0657366) q[3];
sx q[3];
rz(-0.64801884) q[3];
sx q[3];
rz(0.59504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6512063) q[0];
sx q[0];
rz(-0.47845978) q[0];
sx q[0];
rz(-1.867021) q[0];
rz(1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(-1.1638181) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5620016) q[0];
sx q[0];
rz(-2.0193856) q[0];
sx q[0];
rz(-2.6278087) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6357119) q[2];
sx q[2];
rz(-1.9569279) q[2];
sx q[2];
rz(0.94368151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5806942) q[1];
sx q[1];
rz(-0.95798641) q[1];
sx q[1];
rz(-1.3152907) q[1];
rz(-pi) q[2];
rz(-1.9979565) q[3];
sx q[3];
rz(-1.8619976) q[3];
sx q[3];
rz(-1.8607651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.290648) q[2];
sx q[2];
rz(-0.023357563) q[2];
sx q[2];
rz(1.5811496) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-2.0932308) q[3];
sx q[3];
rz(2.7010664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59349638) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(0.12259677) q[0];
rz(-2.3444046) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(3.0534993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6428292) q[0];
sx q[0];
rz(-0.82422148) q[0];
sx q[0];
rz(-0.3513263) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7983059) q[2];
sx q[2];
rz(-1.6517793) q[2];
sx q[2];
rz(-2.0264268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8338203) q[1];
sx q[1];
rz(-1.7227748) q[1];
sx q[1];
rz(0.13843468) q[1];
rz(-pi) q[2];
rz(0.54666211) q[3];
sx q[3];
rz(-0.61184363) q[3];
sx q[3];
rz(0.33239588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(0.49805182) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(-2.3948885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8525456) q[0];
sx q[0];
rz(-0.15327029) q[0];
sx q[0];
rz(2.9943941) q[0];
rz(2.4920801) q[1];
sx q[1];
rz(-1.6308866) q[1];
sx q[1];
rz(1.6726327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2788852) q[0];
sx q[0];
rz(-1.577958) q[0];
sx q[0];
rz(1.347752) q[0];
rz(-pi) q[1];
rz(0.18314731) q[2];
sx q[2];
rz(-1.9948655) q[2];
sx q[2];
rz(1.1893502) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33002871) q[1];
sx q[1];
rz(-2.4061476) q[1];
sx q[1];
rz(1.0388908) q[1];
x q[2];
rz(-2.0135512) q[3];
sx q[3];
rz(-1.0337794) q[3];
sx q[3];
rz(-0.90195105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2088251) q[2];
sx q[2];
rz(-2.6471477) q[2];
sx q[2];
rz(-0.24173582) q[2];
rz(0.73511165) q[3];
sx q[3];
rz(-1.3999516) q[3];
sx q[3];
rz(-0.66489768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46297896) q[0];
sx q[0];
rz(-1.047387) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(2.2562476) q[1];
sx q[1];
rz(-1.0095936) q[1];
sx q[1];
rz(-2.5811894) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40582244) q[0];
sx q[0];
rz(-2.5955831) q[0];
sx q[0];
rz(-1.2237134) q[0];
rz(-pi) q[1];
rz(0.29717314) q[2];
sx q[2];
rz(-1.783833) q[2];
sx q[2];
rz(-1.7293255) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.082277) q[1];
sx q[1];
rz(-0.53566414) q[1];
sx q[1];
rz(-1.0470864) q[1];
x q[2];
rz(-2.3123829) q[3];
sx q[3];
rz(-2.5253339) q[3];
sx q[3];
rz(-0.63502888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6178599) q[2];
sx q[2];
rz(-0.79404074) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(-2.2023885) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(-3.052616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27696779) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(2.4846039) q[0];
rz(-2.9981546) q[1];
sx q[1];
rz(-1.5851494) q[1];
sx q[1];
rz(-2.4030446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26510942) q[0];
sx q[0];
rz(-2.4584157) q[0];
sx q[0];
rz(0.18355145) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0792909) q[2];
sx q[2];
rz(-1.7704983) q[2];
sx q[2];
rz(0.26003597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2354898) q[1];
sx q[1];
rz(-0.74550438) q[1];
sx q[1];
rz(-2.0447313) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62885999) q[3];
sx q[3];
rz(-2.326528) q[3];
sx q[3];
rz(-0.12303837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8831545) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(-3.1385885) q[2];
rz(-2.352412) q[3];
sx q[3];
rz(-0.3843669) q[3];
sx q[3];
rz(1.167231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0810735) q[0];
sx q[0];
rz(-2.7923212) q[0];
sx q[0];
rz(2.9935167) q[0];
rz(-2.608346) q[1];
sx q[1];
rz(-0.24594578) q[1];
sx q[1];
rz(1.6291133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5083744) q[0];
sx q[0];
rz(-2.0323159) q[0];
sx q[0];
rz(1.2715696) q[0];
x q[1];
rz(0.61818168) q[2];
sx q[2];
rz(-0.91518171) q[2];
sx q[2];
rz(-0.55936343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.80485) q[1];
sx q[1];
rz(-2.7798493) q[1];
sx q[1];
rz(1.2233578) q[1];
rz(-pi) q[2];
rz(1.809172) q[3];
sx q[3];
rz(-0.73638232) q[3];
sx q[3];
rz(-3.0958946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6988354) q[2];
sx q[2];
rz(-1.4407225) q[2];
sx q[2];
rz(0.093078144) q[2];
rz(3.0794411) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208948) q[0];
sx q[0];
rz(-2.5466205) q[0];
sx q[0];
rz(2.4769532) q[0];
rz(-1.7918034) q[1];
sx q[1];
rz(-0.87769687) q[1];
sx q[1];
rz(0.68914366) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5377602) q[0];
sx q[0];
rz(-1.6869154) q[0];
sx q[0];
rz(-2.8962747) q[0];
x q[1];
rz(-0.49361009) q[2];
sx q[2];
rz(-0.94075946) q[2];
sx q[2];
rz(0.64865998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.047235977) q[1];
sx q[1];
rz(-0.34730372) q[1];
sx q[1];
rz(-0.44012614) q[1];
x q[2];
rz(1.6313305) q[3];
sx q[3];
rz(-1.7352967) q[3];
sx q[3];
rz(-1.1806844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.050921116) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(1.3760759) q[2];
rz(1.9164267) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(-3.0998949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10251481) q[0];
sx q[0];
rz(-0.044487655) q[0];
sx q[0];
rz(0.59120375) q[0];
rz(-2.8996331) q[1];
sx q[1];
rz(-2.0457025) q[1];
sx q[1];
rz(0.42344365) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4225938) q[0];
sx q[0];
rz(-2.7072577) q[0];
sx q[0];
rz(0.78311626) q[0];
x q[1];
rz(1.6744711) q[2];
sx q[2];
rz(-1.0846234) q[2];
sx q[2];
rz(2.1882265) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4462699) q[1];
sx q[1];
rz(-2.5188418) q[1];
sx q[1];
rz(2.0957309) q[1];
x q[2];
rz(-0.54355346) q[3];
sx q[3];
rz(-2.4494236) q[3];
sx q[3];
rz(-2.8080733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5847136) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(-1.1663743) q[2];
rz(-1.1276468) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(2.6166272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135076) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(0.67310131) q[0];
rz(-0.85064864) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(-0.32352111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7009525) q[0];
sx q[0];
rz(-1.3214759) q[0];
sx q[0];
rz(-1.9560157) q[0];
rz(-pi) q[1];
rz(-2.0490626) q[2];
sx q[2];
rz(-1.8169122) q[2];
sx q[2];
rz(-2.9778632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6767576) q[1];
sx q[1];
rz(-0.82993648) q[1];
sx q[1];
rz(-0.45758943) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1296842) q[3];
sx q[3];
rz(-2.5577196) q[3];
sx q[3];
rz(-1.5990822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(0.072048135) q[2];
rz(-1.0142903) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2608248) q[0];
sx q[0];
rz(-1.7171971) q[0];
sx q[0];
rz(1.1203753) q[0];
rz(-0.33521677) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(-0.14090385) q[2];
sx q[2];
rz(-1.0959411) q[2];
sx q[2];
rz(-0.68077722) q[2];
rz(2.9670197) q[3];
sx q[3];
rz(-2.9963507) q[3];
sx q[3];
rz(0.33490845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
