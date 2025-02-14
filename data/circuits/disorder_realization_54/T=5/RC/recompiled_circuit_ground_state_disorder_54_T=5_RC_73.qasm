OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(3.1867653) q[0];
sx q[0];
rz(10.09633) q[0];
rz(-0.99611941) q[1];
sx q[1];
rz(-2.3528407) q[1];
sx q[1];
rz(2.8532343) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0147284) q[0];
sx q[0];
rz(-1.607195) q[0];
sx q[0];
rz(1.3728736) q[0];
x q[1];
rz(0.66418437) q[2];
sx q[2];
rz(-1.786288) q[2];
sx q[2];
rz(1.2601978) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25569281) q[1];
sx q[1];
rz(-1.608537) q[1];
sx q[1];
rz(2.6934212) q[1];
rz(-pi) q[2];
rz(2.6872271) q[3];
sx q[3];
rz(-1.7011257) q[3];
sx q[3];
rz(-2.5443175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70197376) q[2];
sx q[2];
rz(-0.34586033) q[2];
sx q[2];
rz(0.71887476) q[2];
rz(1.4465205) q[3];
sx q[3];
rz(-1.6547763) q[3];
sx q[3];
rz(-2.1046624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0816536) q[0];
sx q[0];
rz(-0.43123284) q[0];
sx q[0];
rz(-0.28847873) q[0];
rz(0.55229315) q[1];
sx q[1];
rz(-1.0918795) q[1];
sx q[1];
rz(-1.1757895) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4058454) q[0];
sx q[0];
rz(-1.7205392) q[0];
sx q[0];
rz(-1.2558054) q[0];
rz(-pi) q[1];
rz(-0.11643683) q[2];
sx q[2];
rz(-1.7431362) q[2];
sx q[2];
rz(-0.72067537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.9811444) q[1];
sx q[1];
rz(-1.6949883) q[1];
sx q[1];
rz(0.9915413) q[1];
x q[2];
rz(3.1404488) q[3];
sx q[3];
rz(-1.5387156) q[3];
sx q[3];
rz(2.6828121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4484078) q[2];
sx q[2];
rz(-1.8790481) q[2];
sx q[2];
rz(-0.022424879) q[2];
rz(1.4261931) q[3];
sx q[3];
rz(-1.8636999) q[3];
sx q[3];
rz(-0.9001596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65396032) q[0];
sx q[0];
rz(-1.910169) q[0];
sx q[0];
rz(-2.3768429) q[0];
rz(-1.359831) q[1];
sx q[1];
rz(-1.2218852) q[1];
sx q[1];
rz(1.8870032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0645197) q[0];
sx q[0];
rz(-1.700634) q[0];
sx q[0];
rz(-0.10420756) q[0];
rz(-pi) q[1];
rz(-2.8171982) q[2];
sx q[2];
rz(-1.430871) q[2];
sx q[2];
rz(3.0382699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65605983) q[1];
sx q[1];
rz(-2.0295709) q[1];
sx q[1];
rz(-1.2478254) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34864254) q[3];
sx q[3];
rz(-2.7789634) q[3];
sx q[3];
rz(1.8809821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7356871) q[2];
sx q[2];
rz(-2.3730706) q[2];
sx q[2];
rz(-3.1206257) q[2];
rz(1.5812801) q[3];
sx q[3];
rz(-1.0617278) q[3];
sx q[3];
rz(2.235152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2727994) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(-1.3695166) q[0];
rz(2.4319793) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(-2.4443464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6402138) q[0];
sx q[0];
rz(-1.9778628) q[0];
sx q[0];
rz(-1.2329007) q[0];
rz(-pi) q[1];
rz(2.4575649) q[2];
sx q[2];
rz(-2.1993756) q[2];
sx q[2];
rz(2.7492439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5647506) q[1];
sx q[1];
rz(-0.71755845) q[1];
sx q[1];
rz(-2.0400042) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9696021) q[3];
sx q[3];
rz(-1.3376457) q[3];
sx q[3];
rz(-1.5267717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-2.1944025) q[2];
sx q[2];
rz(-2.4626125) q[2];
rz(-1.125157) q[3];
sx q[3];
rz(-1.9966639) q[3];
sx q[3];
rz(-0.2230491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.090002447) q[0];
sx q[0];
rz(-1.9380049) q[0];
sx q[0];
rz(-3.0644655) q[0];
rz(-1.1513101) q[1];
sx q[1];
rz(-1.5733893) q[1];
sx q[1];
rz(0.77879771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6662731) q[0];
sx q[0];
rz(-0.059564807) q[0];
sx q[0];
rz(2.2044529) q[0];
rz(-pi) q[1];
rz(-2.0955032) q[2];
sx q[2];
rz(-2.8150442) q[2];
sx q[2];
rz(-2.7401217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2768138) q[1];
sx q[1];
rz(-1.9052231) q[1];
sx q[1];
rz(0.2590551) q[1];
rz(-pi) q[2];
rz(-3.0833427) q[3];
sx q[3];
rz(-0.16213972) q[3];
sx q[3];
rz(-1.3221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6246346) q[2];
sx q[2];
rz(-1.7004852) q[2];
sx q[2];
rz(-1.1901633) q[2];
rz(-0.55365753) q[3];
sx q[3];
rz(-2.5167969) q[3];
sx q[3];
rz(2.4042118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5364285) q[0];
sx q[0];
rz(-1.9909415) q[0];
sx q[0];
rz(0.27106699) q[0];
rz(1.0003264) q[1];
sx q[1];
rz(-1.9258291) q[1];
sx q[1];
rz(2.2727374) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9199182) q[0];
sx q[0];
rz(-0.98548792) q[0];
sx q[0];
rz(1.0006389) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0539114) q[2];
sx q[2];
rz(-0.73747915) q[2];
sx q[2];
rz(-0.35754851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.083243283) q[1];
sx q[1];
rz(-1.740137) q[1];
sx q[1];
rz(-1.2178628) q[1];
rz(-pi) q[2];
rz(2.6713761) q[3];
sx q[3];
rz(-2.3002671) q[3];
sx q[3];
rz(-2.7470392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.038593682) q[2];
sx q[2];
rz(-0.30240348) q[2];
sx q[2];
rz(-1.137286) q[2];
rz(0.31050995) q[3];
sx q[3];
rz(-0.91028428) q[3];
sx q[3];
rz(0.92946068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0248658) q[0];
sx q[0];
rz(-1.4288582) q[0];
sx q[0];
rz(-1.0799991) q[0];
rz(-2.179821) q[1];
sx q[1];
rz(-0.56517833) q[1];
sx q[1];
rz(0.22294179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7540383) q[0];
sx q[0];
rz(-1.5359274) q[0];
sx q[0];
rz(0.0074444093) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2193394) q[2];
sx q[2];
rz(-1.8195565) q[2];
sx q[2];
rz(-2.3671248) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2126951) q[1];
sx q[1];
rz(-1.2948117) q[1];
sx q[1];
rz(2.3174965) q[1];
x q[2];
rz(-3.1268397) q[3];
sx q[3];
rz(-2.1665769) q[3];
sx q[3];
rz(-2.9500913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96233931) q[2];
sx q[2];
rz(-2.8620359) q[2];
sx q[2];
rz(1.9635828) q[2];
rz(-2.8152605) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(2.0012205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9137664) q[0];
sx q[0];
rz(-0.95781177) q[0];
sx q[0];
rz(-0.80818278) q[0];
rz(-2.783964) q[1];
sx q[1];
rz(-1.6817776) q[1];
sx q[1];
rz(-2.0590032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1515231) q[0];
sx q[0];
rz(-2.4381579) q[0];
sx q[0];
rz(-1.3057054) q[0];
rz(-pi) q[1];
rz(-0.2530667) q[2];
sx q[2];
rz(-2.1089206) q[2];
sx q[2];
rz(0.62935621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1114166) q[1];
sx q[1];
rz(-1.7716265) q[1];
sx q[1];
rz(2.3371731) q[1];
x q[2];
rz(-2.0478422) q[3];
sx q[3];
rz(-1.9821321) q[3];
sx q[3];
rz(0.75440948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.064934405) q[2];
sx q[2];
rz(-0.4883464) q[2];
sx q[2];
rz(-1.9261544) q[2];
rz(-2.2293034) q[3];
sx q[3];
rz(-0.89022294) q[3];
sx q[3];
rz(2.598855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.2927581) q[0];
sx q[0];
rz(-2.500535) q[0];
sx q[0];
rz(0.42375281) q[0];
rz(1.9641701) q[1];
sx q[1];
rz(-0.91554987) q[1];
sx q[1];
rz(0.37568572) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5510124) q[0];
sx q[0];
rz(-2.9504021) q[0];
sx q[0];
rz(-1.952233) q[0];
rz(-pi) q[1];
rz(-2.9027391) q[2];
sx q[2];
rz(-1.544853) q[2];
sx q[2];
rz(-2.8544665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5392307) q[1];
sx q[1];
rz(-1.5054387) q[1];
sx q[1];
rz(2.0728302) q[1];
rz(-pi) q[2];
rz(-0.71640941) q[3];
sx q[3];
rz(-1.3103974) q[3];
sx q[3];
rz(-1.474787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11289135) q[2];
sx q[2];
rz(-1.6501004) q[2];
sx q[2];
rz(-2.4493307) q[2];
rz(-2.4568457) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(-3.0013066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58458644) q[0];
sx q[0];
rz(-2.1742915) q[0];
sx q[0];
rz(-1.517357) q[0];
rz(-1.5380305) q[1];
sx q[1];
rz(-1.7944733) q[1];
sx q[1];
rz(-1.921152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6950078) q[0];
sx q[0];
rz(-1.5112226) q[0];
sx q[0];
rz(-2.1936962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41238971) q[2];
sx q[2];
rz(-0.53942108) q[2];
sx q[2];
rz(-0.88482761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5511849) q[1];
sx q[1];
rz(-0.77279323) q[1];
sx q[1];
rz(2.7332952) q[1];
rz(-0.73699215) q[3];
sx q[3];
rz(-1.7919645) q[3];
sx q[3];
rz(-1.2545409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5662235) q[2];
sx q[2];
rz(-2.1341925) q[2];
sx q[2];
rz(2.2929906) q[2];
rz(-2.1620915) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(-0.037467329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10850567) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(2.5820844) q[1];
sx q[1];
rz(-0.63897501) q[1];
sx q[1];
rz(-0.41313304) q[1];
rz(-1.3146567) q[2];
sx q[2];
rz(-1.8825681) q[2];
sx q[2];
rz(-2.067461) q[2];
rz(0.57942617) q[3];
sx q[3];
rz(-1.5204932) q[3];
sx q[3];
rz(1.6411171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
