OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.644637167453766) q[0];
sx q[0];
rz(3.73124638398225) q[0];
sx q[0];
rz(10.51242015361) q[0];
rz(-0.0154068693518639) q[1];
sx q[1];
rz(3.53137305577333) q[1];
sx q[1];
rz(11.4021472692411) q[1];
cx q[1],q[0];
rz(1.31226599216461) q[0];
sx q[0];
rz(3.2716820110851) q[0];
sx q[0];
rz(9.46123744770094) q[0];
rz(0.340163767337799) q[2];
sx q[2];
rz(5.37469926674897) q[2];
sx q[2];
rz(10.6993901491086) q[2];
cx q[2],q[1];
rz(2.0014181137085) q[1];
sx q[1];
rz(4.92599562008912) q[1];
sx q[1];
rz(9.95533547400638) q[1];
rz(1.4948388338089) q[3];
sx q[3];
rz(4.83004847367341) q[3];
sx q[3];
rz(8.67124412058994) q[3];
cx q[3],q[2];
rz(1.24141335487366) q[2];
sx q[2];
rz(3.70327756007249) q[2];
sx q[2];
rz(10.0708110690038) q[2];
rz(-0.255204051733017) q[3];
sx q[3];
rz(3.59866953094537) q[3];
sx q[3];
rz(10.8193944454114) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.0274885892868) q[0];
sx q[0];
rz(2.80599671800668) q[0];
sx q[0];
rz(10.2317976713102) q[0];
rz(-1.45781600475311) q[1];
sx q[1];
rz(2.5648262818628) q[1];
sx q[1];
rz(11.2225429773252) q[1];
cx q[1],q[0];
rz(1.23490738868713) q[0];
sx q[0];
rz(3.46865648229653) q[0];
sx q[0];
rz(10.5083954095761) q[0];
rz(1.15986669063568) q[2];
sx q[2];
rz(4.183537991839) q[2];
sx q[2];
rz(7.91131732463046) q[2];
cx q[2],q[1];
rz(2.46494841575623) q[1];
sx q[1];
rz(2.05687716801698) q[1];
sx q[1];
rz(10.1193327069203) q[1];
rz(1.61557722091675) q[3];
sx q[3];
rz(3.98402258952195) q[3];
sx q[3];
rz(10.8090446948926) q[3];
cx q[3],q[2];
rz(0.96062159538269) q[2];
sx q[2];
rz(4.21169939835603) q[2];
sx q[2];
rz(9.23105118273898) q[2];
rz(1.52420294284821) q[3];
sx q[3];
rz(2.38348552783067) q[3];
sx q[3];
rz(9.60352856516048) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.52095651626587) q[0];
sx q[0];
rz(2.90363235970075) q[0];
sx q[0];
rz(10.4008798956792) q[0];
rz(2.25282645225525) q[1];
sx q[1];
rz(2.39753940899903) q[1];
sx q[1];
rz(9.25173065661594) q[1];
cx q[1],q[0];
rz(0.425874680280685) q[0];
sx q[0];
rz(3.10396553774411) q[0];
sx q[0];
rz(9.3499938532631) q[0];
rz(0.125421807169914) q[2];
sx q[2];
rz(4.19516983826692) q[2];
sx q[2];
rz(9.77461708187267) q[2];
cx q[2],q[1];
rz(1.08066070079803) q[1];
sx q[1];
rz(2.0736040194803) q[1];
sx q[1];
rz(9.23541869818374) q[1];
rz(1.2552433013916) q[3];
sx q[3];
rz(3.47951919038827) q[3];
sx q[3];
rz(9.03779960273906) q[3];
cx q[3],q[2];
rz(-1.13910698890686) q[2];
sx q[2];
rz(1.3287012894922) q[2];
sx q[2];
rz(11.852496123306) q[2];
rz(-0.21126963198185) q[3];
sx q[3];
rz(4.17887643178041) q[3];
sx q[3];
rz(8.85941407679721) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.367703765630722) q[0];
sx q[0];
rz(4.32508209546144) q[0];
sx q[0];
rz(10.55756685733) q[0];
rz(-0.471381306648254) q[1];
sx q[1];
rz(4.43551305134828) q[1];
sx q[1];
rz(8.99376440643474) q[1];
cx q[1],q[0];
rz(0.339922666549683) q[0];
sx q[0];
rz(3.3819801231199) q[0];
sx q[0];
rz(10.5024034738462) q[0];
rz(0.75881826877594) q[2];
sx q[2];
rz(5.78013339837129) q[2];
sx q[2];
rz(8.69630709885761) q[2];
cx q[2],q[1];
rz(1.18725943565369) q[1];
sx q[1];
rz(2.10128393967683) q[1];
sx q[1];
rz(9.775710707895) q[1];
rz(-0.621638119220734) q[3];
sx q[3];
rz(4.49007955391938) q[3];
sx q[3];
rz(10.1779879689138) q[3];
cx q[3],q[2];
rz(1.62125778198242) q[2];
sx q[2];
rz(3.88112756808335) q[2];
sx q[2];
rz(9.83106035589381) q[2];
rz(1.20643651485443) q[3];
sx q[3];
rz(2.04835692246492) q[3];
sx q[3];
rz(11.9796976804654) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.632211446762085) q[0];
sx q[0];
rz(4.01726207335527) q[0];
sx q[0];
rz(9.75111523865863) q[0];
rz(2.18747663497925) q[1];
sx q[1];
rz(2.49311444361741) q[1];
sx q[1];
rz(9.40125315486594) q[1];
cx q[1],q[0];
rz(-0.017880579456687) q[0];
sx q[0];
rz(2.57662257750566) q[0];
sx q[0];
rz(10.1739573240201) q[0];
rz(2.82511878013611) q[2];
sx q[2];
rz(4.07883325417573) q[2];
sx q[2];
rz(8.75180194377109) q[2];
cx q[2],q[1];
rz(0.702908039093018) q[1];
sx q[1];
rz(3.59506181080873) q[1];
sx q[1];
rz(10.1221158265988) q[1];
rz(1.13773167133331) q[3];
sx q[3];
rz(3.46522921522195) q[3];
sx q[3];
rz(8.10466990470096) q[3];
cx q[3],q[2];
rz(0.247133687138557) q[2];
sx q[2];
rz(3.60025406082208) q[2];
sx q[2];
rz(10.0886952042501) q[2];
rz(-0.890200555324554) q[3];
sx q[3];
rz(4.69084134896333) q[3];
sx q[3];
rz(9.43748336787477) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.35060548782349) q[0];
sx q[0];
rz(3.57302144368226) q[0];
sx q[0];
rz(9.75560129284068) q[0];
rz(0.362589687108994) q[1];
sx q[1];
rz(4.32090106804902) q[1];
sx q[1];
rz(9.94052538870975) q[1];
cx q[1],q[0];
rz(-0.121304854750633) q[0];
sx q[0];
rz(3.46563220222528) q[0];
sx q[0];
rz(8.68552759884998) q[0];
rz(-0.976769685745239) q[2];
sx q[2];
rz(4.55451742013032) q[2];
sx q[2];
rz(11.2342118978421) q[2];
cx q[2],q[1];
rz(-0.679298996925354) q[1];
sx q[1];
rz(4.0223418196016) q[1];
sx q[1];
rz(11.0959665536801) q[1];
rz(-0.132388740777969) q[3];
sx q[3];
rz(3.92087206442887) q[3];
sx q[3];
rz(9.17326763867542) q[3];
cx q[3],q[2];
rz(0.350468784570694) q[2];
sx q[2];
rz(1.45012691815431) q[2];
sx q[2];
rz(8.61725852488681) q[2];
rz(0.098585419356823) q[3];
sx q[3];
rz(5.38512221177156) q[3];
sx q[3];
rz(10.5175883531491) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.438078552484512) q[0];
sx q[0];
rz(3.71483156283433) q[0];
sx q[0];
rz(9.87430760859653) q[0];
rz(2.45641565322876) q[1];
sx q[1];
rz(3.54919842083985) q[1];
sx q[1];
rz(10.0442464709203) q[1];
cx q[1],q[0];
rz(2.02800583839417) q[0];
sx q[0];
rz(5.02186194260652) q[0];
sx q[0];
rz(9.45449659823581) q[0];
rz(1.04256045818329) q[2];
sx q[2];
rz(1.52914574940736) q[2];
sx q[2];
rz(10.741853094093) q[2];
cx q[2],q[1];
rz(0.0723143145442009) q[1];
sx q[1];
rz(5.3855072577768) q[1];
sx q[1];
rz(9.85901997088596) q[1];
rz(0.731605589389801) q[3];
sx q[3];
rz(3.41424179275567) q[3];
sx q[3];
rz(10.38817689418) q[3];
cx q[3],q[2];
rz(0.56457906961441) q[2];
sx q[2];
rz(4.83422199090058) q[2];
sx q[2];
rz(9.22676866351768) q[2];
rz(1.36245822906494) q[3];
sx q[3];
rz(2.84153157671029) q[3];
sx q[3];
rz(8.72438189982578) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.144415974617) q[0];
sx q[0];
rz(3.76348284085328) q[0];
sx q[0];
rz(11.280811405174) q[0];
rz(-0.276649117469788) q[1];
sx q[1];
rz(4.51028338273103) q[1];
sx q[1];
rz(9.31616266667053) q[1];
cx q[1],q[0];
rz(0.40688693523407) q[0];
sx q[0];
rz(2.63317749102647) q[0];
sx q[0];
rz(9.5935170262973) q[0];
rz(-0.109933964908123) q[2];
sx q[2];
rz(3.76846626599366) q[2];
sx q[2];
rz(9.3813173904936) q[2];
cx q[2],q[1];
rz(1.36563551425934) q[1];
sx q[1];
rz(4.49808338482911) q[1];
sx q[1];
rz(10.1448911189954) q[1];
rz(0.386445552110672) q[3];
sx q[3];
rz(4.59579375584657) q[3];
sx q[3];
rz(10.4618609905164) q[3];
cx q[3],q[2];
rz(0.308399766683578) q[2];
sx q[2];
rz(4.33822110493714) q[2];
sx q[2];
rz(11.2921151876371) q[2];
rz(1.3537003993988) q[3];
sx q[3];
rz(3.58089852531488) q[3];
sx q[3];
rz(8.5581943154256) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.152477905154228) q[0];
sx q[0];
rz(3.86725226243074) q[0];
sx q[0];
rz(9.3042231336157) q[0];
rz(-0.741479158401489) q[1];
sx q[1];
rz(4.34560671647126) q[1];
sx q[1];
rz(9.50046484022542) q[1];
cx q[1],q[0];
rz(0.736294329166412) q[0];
sx q[0];
rz(4.40672305424745) q[0];
sx q[0];
rz(9.59690608679458) q[0];
rz(0.15472137928009) q[2];
sx q[2];
rz(3.60795632203157) q[2];
sx q[2];
rz(11.0613744020383) q[2];
cx q[2],q[1];
rz(-0.329667419195175) q[1];
sx q[1];
rz(3.79823050101335) q[1];
sx q[1];
rz(11.0369765520017) q[1];
rz(-1.30660605430603) q[3];
sx q[3];
rz(3.75909987290437) q[3];
sx q[3];
rz(10.1914726853292) q[3];
cx q[3],q[2];
rz(-0.0724482983350754) q[2];
sx q[2];
rz(4.02439502080018) q[2];
sx q[2];
rz(9.77111182212039) q[2];
rz(0.945414662361145) q[3];
sx q[3];
rz(4.46411481698091) q[3];
sx q[3];
rz(10.3719277739446) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.27461767196655) q[0];
sx q[0];
rz(4.3397468646341) q[0];
sx q[0];
rz(9.87266681193515) q[0];
rz(0.92944324016571) q[1];
sx q[1];
rz(4.54106167157228) q[1];
sx q[1];
rz(8.75026050805255) q[1];
cx q[1],q[0];
rz(0.697929918766022) q[0];
sx q[0];
rz(3.95723954041535) q[0];
sx q[0];
rz(9.70112306474849) q[0];
rz(0.506979942321777) q[2];
sx q[2];
rz(3.35339032311971) q[2];
sx q[2];
rz(10.7422660350721) q[2];
cx q[2],q[1];
rz(1.06219398975372) q[1];
sx q[1];
rz(3.31498950918252) q[1];
sx q[1];
rz(9.4313231007117) q[1];
rz(-0.418728440999985) q[3];
sx q[3];
rz(4.6434567292505) q[3];
sx q[3];
rz(10.9381489515226) q[3];
cx q[3],q[2];
rz(0.210053607821465) q[2];
sx q[2];
rz(4.35117462475831) q[2];
sx q[2];
rz(9.69217372535869) q[2];
rz(1.10726547241211) q[3];
sx q[3];
rz(3.92406991322572) q[3];
sx q[3];
rz(10.2737697720449) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.258765906095505) q[0];
sx q[0];
rz(3.6280629654699) q[0];
sx q[0];
rz(9.31369948982402) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-2.52006006240845) q[1];
sx q[1];
rz(3.66290161212022) q[1];
sx q[1];
rz(12.7334804296414) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.152059853076935) q[2];
sx q[2];
rz(4.57758048375184) q[2];
sx q[2];
rz(10.9032832145612) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.02916991710663) q[3];
sx q[3];
rz(3.62468081911141) q[3];
sx q[3];
rz(9.28242153524562) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
