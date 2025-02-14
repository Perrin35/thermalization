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
rz(2.661854) q[0];
sx q[0];
rz(5.148611) q[0];
sx q[0];
rz(10.892286) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(-2.1958513) q[1];
sx q[1];
rz(-2.6374964) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5382982) q[0];
sx q[0];
rz(-2.3497938) q[0];
sx q[0];
rz(-2.3843423) q[0];
rz(0.43457793) q[2];
sx q[2];
rz(-0.38056669) q[2];
sx q[2];
rz(-1.6004882) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13611804) q[1];
sx q[1];
rz(-2.1393993) q[1];
sx q[1];
rz(-0.75638812) q[1];
rz(0.33511843) q[3];
sx q[3];
rz(-1.5247282) q[3];
sx q[3];
rz(0.81460458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(-0.18623713) q[2];
rz(1.7765744) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(1.9369102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-0.052074281) q[0];
rz(2.1049818) q[1];
sx q[1];
rz(-0.60984817) q[1];
sx q[1];
rz(0.61526543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0947845) q[0];
sx q[0];
rz(-2.2333217) q[0];
sx q[0];
rz(-2.725968) q[0];
rz(-pi) q[1];
rz(-2.6411166) q[2];
sx q[2];
rz(-2.8722974) q[2];
sx q[2];
rz(0.057502086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2959472) q[1];
sx q[1];
rz(-1.6270279) q[1];
sx q[1];
rz(-0.91744951) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7795981) q[3];
sx q[3];
rz(-2.3376645) q[3];
sx q[3];
rz(1.7957791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.44876862) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(1.2358865) q[2];
rz(1.1290733) q[3];
sx q[3];
rz(-2.2344507) q[3];
sx q[3];
rz(-0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579984) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(0.54028571) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(1.6859863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42403635) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(-2.9633029) q[0];
x q[1];
rz(-1.6605145) q[2];
sx q[2];
rz(-1.3680581) q[2];
sx q[2];
rz(-0.81344542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8639455) q[1];
sx q[1];
rz(-1.2182115) q[1];
sx q[1];
rz(-0.90971281) q[1];
rz(2.0469401) q[3];
sx q[3];
rz(-0.75679251) q[3];
sx q[3];
rz(-1.4534637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-0.33438412) q[2];
sx q[2];
rz(1.4462659) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-0.97508591) q[3];
sx q[3];
rz(1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0541662) q[0];
sx q[0];
rz(-0.26987258) q[0];
sx q[0];
rz(0.42359459) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(0.097600309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256297) q[0];
sx q[0];
rz(-2.0920252) q[0];
sx q[0];
rz(2.002191) q[0];
rz(-2.9984442) q[2];
sx q[2];
rz(-2.2586056) q[2];
sx q[2];
rz(-2.0542415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.017189055) q[1];
sx q[1];
rz(-1.9162971) q[1];
sx q[1];
rz(-0.33331897) q[1];
rz(-pi) q[2];
rz(1.9097435) q[3];
sx q[3];
rz(-2.2710861) q[3];
sx q[3];
rz(0.46745121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5264954) q[2];
sx q[2];
rz(-2.6690833) q[2];
sx q[2];
rz(0.16461593) q[2];
rz(0.047867157) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(-0.78711787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98016244) q[0];
sx q[0];
rz(-2.6949368) q[0];
sx q[0];
rz(-1.5572146) q[0];
rz(0.36852512) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(-1.535086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6045348) q[0];
sx q[0];
rz(-0.7365444) q[0];
sx q[0];
rz(-2.3307072) q[0];
x q[1];
rz(-2.317165) q[2];
sx q[2];
rz(-1.8608837) q[2];
sx q[2];
rz(0.96765358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45101825) q[1];
sx q[1];
rz(-1.3214045) q[1];
sx q[1];
rz(-1.0728157) q[1];
rz(-0.078952958) q[3];
sx q[3];
rz(-1.8584849) q[3];
sx q[3];
rz(-2.1514078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(2.0965915) q[2];
rz(0.16799489) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(-2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217693) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(-1.2044915) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-0.81356994) q[1];
sx q[1];
rz(2.4258851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2432118) q[0];
sx q[0];
rz(-2.1247433) q[0];
sx q[0];
rz(2.5701163) q[0];
rz(-pi) q[1];
rz(-0.44797051) q[2];
sx q[2];
rz(-1.391419) q[2];
sx q[2];
rz(-0.57283869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5454272) q[1];
sx q[1];
rz(-1.6837032) q[1];
sx q[1];
rz(2.4775056) q[1];
rz(-pi) q[2];
rz(1.1119858) q[3];
sx q[3];
rz(-1.8355287) q[3];
sx q[3];
rz(-2.2609425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(2.0984207) q[2];
rz(-0.10797524) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(2.030355) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1825948) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(1.2249655) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-0.8871612) q[1];
sx q[1];
rz(1.8909594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179493) q[0];
sx q[0];
rz(-2.3139489) q[0];
sx q[0];
rz(-1.4046937) q[0];
rz(-pi) q[1];
rz(-3.1194889) q[2];
sx q[2];
rz(-1.7689509) q[2];
sx q[2];
rz(-1.0304067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1742184) q[1];
sx q[1];
rz(-1.3506445) q[1];
sx q[1];
rz(-2.8323814) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5851303) q[3];
sx q[3];
rz(-1.5598522) q[3];
sx q[3];
rz(1.3370322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(2.4459548) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(-2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19145963) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(-2.4323442) q[0];
rz(1.5059772) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(1.0848612) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1902311) q[0];
sx q[0];
rz(-1.5133281) q[0];
sx q[0];
rz(-0.60110737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80471595) q[2];
sx q[2];
rz(-1.448481) q[2];
sx q[2];
rz(1.027077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49902713) q[1];
sx q[1];
rz(-2.0303747) q[1];
sx q[1];
rz(-1.2867371) q[1];
rz(-pi) q[2];
rz(-0.011452518) q[3];
sx q[3];
rz(-1.6310777) q[3];
sx q[3];
rz(1.1520916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52013493) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(1.8247068) q[2];
rz(2.3560611) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21996466) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(-0.40837902) q[0];
rz(2.1829103) q[1];
sx q[1];
rz(-2.8125693) q[1];
sx q[1];
rz(-1.5214517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565529) q[0];
sx q[0];
rz(-1.3360177) q[0];
sx q[0];
rz(1.9454369) q[0];
rz(-3.0651363) q[2];
sx q[2];
rz(-2.5707013) q[2];
sx q[2];
rz(0.61797188) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16220763) q[1];
sx q[1];
rz(-1.1623135) q[1];
sx q[1];
rz(-0.86345203) q[1];
rz(0.19100325) q[3];
sx q[3];
rz(-1.2634209) q[3];
sx q[3];
rz(1.5366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7659605) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(2.5755303) q[2];
rz(-1.3426956) q[3];
sx q[3];
rz(-0.415396) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.600243) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(-1.4254697) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(-0.38965449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80042875) q[0];
sx q[0];
rz(-2.6376056) q[0];
sx q[0];
rz(-2.1227073) q[0];
x q[1];
rz(-0.25954982) q[2];
sx q[2];
rz(-2.1786111) q[2];
sx q[2];
rz(1.1961301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13481465) q[1];
sx q[1];
rz(-1.14535) q[1];
sx q[1];
rz(-0.025789217) q[1];
rz(-pi) q[2];
rz(1.2663307) q[3];
sx q[3];
rz(-1.0926477) q[3];
sx q[3];
rz(-1.4451417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9044372) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(-1.944444) q[2];
rz(0.14687471) q[3];
sx q[3];
rz(-0.67320383) q[3];
sx q[3];
rz(0.98744121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.7214397) q[0];
sx q[0];
rz(-1.0316105) q[0];
sx q[0];
rz(1.6461865) q[0];
rz(2.738476) q[1];
sx q[1];
rz(-1.3251726) q[1];
sx q[1];
rz(-1.6641738) q[1];
rz(2.9523136) q[2];
sx q[2];
rz(-0.58313417) q[2];
sx q[2];
rz(-1.492576) q[2];
rz(0.38539376) q[3];
sx q[3];
rz(-1.3968395) q[3];
sx q[3];
rz(-2.2341961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
