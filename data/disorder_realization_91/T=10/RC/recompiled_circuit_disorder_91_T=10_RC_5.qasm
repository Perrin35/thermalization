OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(5.7167238) q[0];
sx q[0];
rz(9.5958435) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1216461) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(-1.6365504) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1034307) q[2];
sx q[2];
rz(-0.51617235) q[2];
sx q[2];
rz(-1.49522) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6273856) q[1];
sx q[1];
rz(-1.5878907) q[1];
sx q[1];
rz(1.3196726) q[1];
x q[2];
rz(0.29403789) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.99130327) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(1.4527028) q[0];
rz(0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(0.30028775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072104134) q[0];
sx q[0];
rz(-0.61394962) q[0];
sx q[0];
rz(0.21013649) q[0];
x q[1];
rz(-0.57057256) q[2];
sx q[2];
rz(-1.8405481) q[2];
sx q[2];
rz(2.3938993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.7784361) q[1];
sx q[1];
rz(-1.8646212) q[1];
sx q[1];
rz(-0.2328965) q[1];
rz(-pi) q[2];
rz(-2.977936) q[3];
sx q[3];
rz(-1.3842584) q[3];
sx q[3];
rz(1.877762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(-2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5791941) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(2.753479) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(0.31633502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280699) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(1.3257922) q[0];
rz(-1.9947617) q[2];
sx q[2];
rz(-2.1018873) q[2];
sx q[2];
rz(2.6742427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0132732) q[1];
sx q[1];
rz(-2.0187906) q[1];
sx q[1];
rz(-2.0112579) q[1];
x q[2];
rz(2.8510677) q[3];
sx q[3];
rz(-1.8968664) q[3];
sx q[3];
rz(0.80929221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0324273) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(-1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.9008824) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-0.8316935) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261525) q[0];
sx q[0];
rz(-0.8988131) q[0];
sx q[0];
rz(-1.6295022) q[0];
x q[1];
rz(-2.1349147) q[2];
sx q[2];
rz(-1.5022105) q[2];
sx q[2];
rz(-2.2401631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6935389) q[1];
sx q[1];
rz(-1.5539907) q[1];
sx q[1];
rz(-1.1225213) q[1];
rz(-pi) q[2];
rz(1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(0.26879877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(-1.303724) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542959) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(3.0412578) q[0];
x q[1];
rz(0.15437834) q[2];
sx q[2];
rz(-1.2882243) q[2];
sx q[2];
rz(-2.4992361) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3225973) q[1];
sx q[1];
rz(-1.4154735) q[1];
sx q[1];
rz(1.1197234) q[1];
x q[2];
rz(-0.13637654) q[3];
sx q[3];
rz(-2.4032421) q[3];
sx q[3];
rz(2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6872528) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(2.5435737) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(-2.0050744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-0.20198527) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.252713) q[0];
sx q[0];
rz(-1.3759334) q[0];
sx q[0];
rz(-1.7478419) q[0];
rz(-pi) q[1];
rz(-1.8385356) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(-2.7051085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0373842) q[1];
sx q[1];
rz(-1.6777778) q[1];
sx q[1];
rz(-1.2030829) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5842651) q[3];
sx q[3];
rz(-1.8263655) q[3];
sx q[3];
rz(-0.42423466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.3827682) q[2];
rz(2.8921195) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(-1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7145342) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(-2.2447341) q[0];
rz(2.8915021) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8539124) q[0];
sx q[0];
rz(-2.2932862) q[0];
sx q[0];
rz(1.393909) q[0];
rz(-1.8115225) q[2];
sx q[2];
rz(-2.4154818) q[2];
sx q[2];
rz(-1.7118529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0301789) q[1];
sx q[1];
rz(-1.3707268) q[1];
sx q[1];
rz(2.5740037) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2956804) q[3];
sx q[3];
rz(-1.0529622) q[3];
sx q[3];
rz(-0.3412316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0573132) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-0.21952595) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.052208386) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(0.6589748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7227104) q[0];
sx q[0];
rz(-1.5960777) q[0];
sx q[0];
rz(-0.77194571) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72139229) q[2];
sx q[2];
rz(-2.112769) q[2];
sx q[2];
rz(-1.1162356) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36775667) q[1];
sx q[1];
rz(-0.78513297) q[1];
sx q[1];
rz(-2.8326041) q[1];
x q[2];
rz(0.068816618) q[3];
sx q[3];
rz(-0.92888442) q[3];
sx q[3];
rz(1.2115692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3808688) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.9314996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8378728) q[0];
sx q[0];
rz(-1.923773) q[0];
sx q[0];
rz(2.9205802) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96610539) q[2];
sx q[2];
rz(-1.3622869) q[2];
sx q[2];
rz(1.0143806) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.853211) q[1];
sx q[1];
rz(-2.2182811) q[1];
sx q[1];
rz(2.7700469) q[1];
rz(2.5217767) q[3];
sx q[3];
rz(-1.8780939) q[3];
sx q[3];
rz(0.67051552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(-0.56662095) q[2];
rz(2.2120655) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(-1.43191) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0678588) q[0];
sx q[0];
rz(-1.6556544) q[0];
sx q[0];
rz(-1.9195062) q[0];
rz(-pi) q[1];
rz(-1.6522371) q[2];
sx q[2];
rz(-1.7310206) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.666115) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(2.75027) q[1];
x q[2];
rz(2.6639789) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(-1.175566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8563103) q[2];
sx q[2];
rz(-0.44390634) q[2];
sx q[2];
rz(2.4463859) q[2];
rz(2.5837512) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0859062) q[0];
sx q[0];
rz(-1.1924556) q[0];
sx q[0];
rz(-2.6299155) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(1.5716626) q[2];
sx q[2];
rz(-1.9067087) q[2];
sx q[2];
rz(-2.9839347) q[2];
rz(-2.3435076) q[3];
sx q[3];
rz(-1.4555664) q[3];
sx q[3];
rz(1.2496787) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
