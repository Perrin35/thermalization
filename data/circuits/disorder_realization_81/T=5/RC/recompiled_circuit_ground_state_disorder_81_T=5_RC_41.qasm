OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(-1.4338926) q[1];
sx q[1];
rz(1.5418928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34719742) q[0];
sx q[0];
rz(-1.3566249) q[0];
sx q[0];
rz(-1.2783575) q[0];
rz(-pi) q[1];
rz(-2.9028106) q[2];
sx q[2];
rz(-0.14792837) q[2];
sx q[2];
rz(1.5419568) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1356537) q[1];
sx q[1];
rz(-2.7126813) q[1];
sx q[1];
rz(-0.072838293) q[1];
x q[2];
rz(2.466731) q[3];
sx q[3];
rz(-2.6913342) q[3];
sx q[3];
rz(-0.78801149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8371007) q[2];
sx q[2];
rz(-1.203275) q[2];
sx q[2];
rz(0.692918) q[2];
rz(-0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(-2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8189341) q[0];
sx q[0];
rz(-2.942473) q[0];
sx q[0];
rz(-2.0767427) q[0];
rz(0.62243593) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(-2.650824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7642794) q[0];
sx q[0];
rz(-1.547033) q[0];
sx q[0];
rz(1.5800493) q[0];
x q[1];
rz(-1.7945917) q[2];
sx q[2];
rz(-1.9836805) q[2];
sx q[2];
rz(-0.041963581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0239531) q[1];
sx q[1];
rz(-2.3815386) q[1];
sx q[1];
rz(-3.0031086) q[1];
rz(-pi) q[2];
rz(1.4014259) q[3];
sx q[3];
rz(-1.6182096) q[3];
sx q[3];
rz(0.58978888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48851594) q[2];
sx q[2];
rz(-1.642903) q[2];
sx q[2];
rz(-2.901279) q[2];
rz(2.6416685) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8847454) q[0];
sx q[0];
rz(-0.95504967) q[0];
sx q[0];
rz(0.98168674) q[0];
rz(-0.49199545) q[1];
sx q[1];
rz(-1.0849846) q[1];
sx q[1];
rz(-1.5031776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5444191) q[0];
sx q[0];
rz(-1.0962631) q[0];
sx q[0];
rz(-2.3719792) q[0];
x q[1];
rz(-1.6257587) q[2];
sx q[2];
rz(-1.9703437) q[2];
sx q[2];
rz(1.1732701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7780077) q[1];
sx q[1];
rz(-2.6202046) q[1];
sx q[1];
rz(-0.85659493) q[1];
rz(-pi) q[2];
rz(-1.7812626) q[3];
sx q[3];
rz(-0.099848824) q[3];
sx q[3];
rz(2.9304402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65695277) q[2];
sx q[2];
rz(-0.52565614) q[2];
sx q[2];
rz(-3.0204115) q[2];
rz(-1.0080522) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693531) q[0];
sx q[0];
rz(-3.1407052) q[0];
sx q[0];
rz(2.6278611) q[0];
rz(0.95798245) q[1];
sx q[1];
rz(-1.1424516) q[1];
sx q[1];
rz(1.4428008) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0009036) q[0];
sx q[0];
rz(-1.1620518) q[0];
sx q[0];
rz(-1.391414) q[0];
rz(-pi) q[1];
rz(1.3495693) q[2];
sx q[2];
rz(-2.0906679) q[2];
sx q[2];
rz(1.9398361) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5123547) q[1];
sx q[1];
rz(-1.6529473) q[1];
sx q[1];
rz(-0.95463971) q[1];
rz(-pi) q[2];
rz(0.12421457) q[3];
sx q[3];
rz(-0.66541568) q[3];
sx q[3];
rz(-1.6197325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55502597) q[2];
sx q[2];
rz(-0.96379605) q[2];
sx q[2];
rz(-1.2653992) q[2];
rz(1.966656) q[3];
sx q[3];
rz(-2.411071) q[3];
sx q[3];
rz(-1.9780212) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555399) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(1.3245921) q[0];
rz(0.96877226) q[1];
sx q[1];
rz(-1.4854393) q[1];
sx q[1];
rz(2.103215) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85920716) q[0];
sx q[0];
rz(-2.7250054) q[0];
sx q[0];
rz(0.4075012) q[0];
rz(-pi) q[1];
rz(2.3275787) q[2];
sx q[2];
rz(-2.2377439) q[2];
sx q[2];
rz(2.2499354) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5165625) q[1];
sx q[1];
rz(-1.6119526) q[1];
sx q[1];
rz(2.7709318) q[1];
x q[2];
rz(2.4828827) q[3];
sx q[3];
rz(-1.1532056) q[3];
sx q[3];
rz(-2.5233248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5021299) q[2];
sx q[2];
rz(-2.837193) q[2];
sx q[2];
rz(-0.89140618) q[2];
rz(-1.8799479) q[3];
sx q[3];
rz(-1.6311878) q[3];
sx q[3];
rz(0.6663028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0618184) q[0];
sx q[0];
rz(-2.6462055) q[0];
sx q[0];
rz(-0.99639446) q[0];
rz(-2.2415316) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(0.17734227) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9422007) q[0];
sx q[0];
rz(-0.15155242) q[0];
sx q[0];
rz(-1.5723537) q[0];
x q[1];
rz(2.2093973) q[2];
sx q[2];
rz(-1.5538486) q[2];
sx q[2];
rz(-1.0826966) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85834938) q[1];
sx q[1];
rz(-0.56811404) q[1];
sx q[1];
rz(-1.0198221) q[1];
x q[2];
rz(-3.0445339) q[3];
sx q[3];
rz(-2.0776111) q[3];
sx q[3];
rz(-0.34891303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.29048723) q[2];
sx q[2];
rz(-1.9600211) q[2];
sx q[2];
rz(-0.30430749) q[2];
rz(-2.815222) q[3];
sx q[3];
rz(-2.5984952) q[3];
sx q[3];
rz(-0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0684763) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(-2.2090744) q[0];
rz(-0.069123507) q[1];
sx q[1];
rz(-2.9239475) q[1];
sx q[1];
rz(-2.1014012) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016259232) q[0];
sx q[0];
rz(-2.3361492) q[0];
sx q[0];
rz(-1.6065434) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97904737) q[2];
sx q[2];
rz(-2.8354916) q[2];
sx q[2];
rz(3.0476168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2582764) q[1];
sx q[1];
rz(-0.81352106) q[1];
sx q[1];
rz(-1.3431647) q[1];
rz(-pi) q[2];
rz(1.0024985) q[3];
sx q[3];
rz(-2.1868984) q[3];
sx q[3];
rz(1.440986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13504623) q[2];
sx q[2];
rz(-0.73837215) q[2];
sx q[2];
rz(-0.73299232) q[2];
rz(2.0586355) q[3];
sx q[3];
rz(-2.0854009) q[3];
sx q[3];
rz(2.1347031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873544) q[0];
sx q[0];
rz(-3.0576958) q[0];
sx q[0];
rz(2.6485637) q[0];
rz(-0.21656491) q[1];
sx q[1];
rz(-1.1410057) q[1];
sx q[1];
rz(-2.1176178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7465377) q[0];
sx q[0];
rz(-0.97486541) q[0];
sx q[0];
rz(-1.5563957) q[0];
rz(1.4570974) q[2];
sx q[2];
rz(-1.2405044) q[2];
sx q[2];
rz(1.8735261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1223381) q[1];
sx q[1];
rz(-0.38568364) q[1];
sx q[1];
rz(-2.7606439) q[1];
rz(-pi) q[2];
rz(-0.29924691) q[3];
sx q[3];
rz(-0.95029921) q[3];
sx q[3];
rz(2.3503691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0852802) q[2];
sx q[2];
rz(-1.6322501) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(0.41302776) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(3.107048) q[0];
rz(3.0870364) q[1];
sx q[1];
rz(-1.1628393) q[1];
sx q[1];
rz(2.3202855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019857865) q[0];
sx q[0];
rz(-2.1424751) q[0];
sx q[0];
rz(1.4522533) q[0];
x q[1];
rz(3.0376833) q[2];
sx q[2];
rz(-0.089251554) q[2];
sx q[2];
rz(-0.98429843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2945613) q[1];
sx q[1];
rz(-0.18974133) q[1];
sx q[1];
rz(-2.4824597) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2230439) q[3];
sx q[3];
rz(-1.9652742) q[3];
sx q[3];
rz(-0.6834417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7811232) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(-0.13295573) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.6229595) q[3];
sx q[3];
rz(2.0558004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0947615) q[0];
sx q[0];
rz(-2.5346041) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(-1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(0.81319317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614118) q[0];
sx q[0];
rz(-1.6488355) q[0];
sx q[0];
rz(-1.5355996) q[0];
rz(-pi) q[1];
rz(0.1172448) q[2];
sx q[2];
rz(-2.239253) q[2];
sx q[2];
rz(0.14456597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.606914) q[1];
sx q[1];
rz(-0.77086222) q[1];
sx q[1];
rz(-2.4669564) q[1];
x q[2];
rz(0.45222262) q[3];
sx q[3];
rz(-1.4503696) q[3];
sx q[3];
rz(2.6129006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64858156) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(-1.0515593) q[2];
rz(-2.1227396) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(-0.65783182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(1.3195994) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(-0.32991275) q[2];
sx q[2];
rz(-2.920426) q[2];
sx q[2];
rz(-2.5830808) q[2];
rz(0.29215688) q[3];
sx q[3];
rz(-2.9360129) q[3];
sx q[3];
rz(2.871411) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
