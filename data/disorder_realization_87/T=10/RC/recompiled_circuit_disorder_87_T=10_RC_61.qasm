OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(0.2015764) q[0];
rz(0.4959313) q[1];
sx q[1];
rz(-0.5402686) q[1];
sx q[1];
rz(2.2044866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21298458) q[0];
sx q[0];
rz(-2.1323418) q[0];
sx q[0];
rz(-2.0748078) q[0];
rz(-pi) q[1];
rz(2.4891698) q[2];
sx q[2];
rz(-0.7235652) q[2];
sx q[2];
rz(-3.0384118) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.571849) q[1];
sx q[1];
rz(-1.7535926) q[1];
sx q[1];
rz(0.079449541) q[1];
rz(3.044371) q[3];
sx q[3];
rz(-1.3029649) q[3];
sx q[3];
rz(0.83084805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9822838) q[2];
sx q[2];
rz(-0.26580492) q[2];
sx q[2];
rz(-3.0554331) q[2];
rz(-0.75749767) q[3];
sx q[3];
rz(-0.75452724) q[3];
sx q[3];
rz(2.0479726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5646097) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(-1.3264867) q[0];
rz(1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87646987) q[0];
sx q[0];
rz(-2.070283) q[0];
sx q[0];
rz(-0.45544099) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6255252) q[2];
sx q[2];
rz(-1.3000254) q[2];
sx q[2];
rz(-2.5004636) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2084864) q[1];
sx q[1];
rz(-2.3098364) q[1];
sx q[1];
rz(-2.0887124) q[1];
rz(-0.80941697) q[3];
sx q[3];
rz(-1.8339694) q[3];
sx q[3];
rz(-2.0078878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0469971) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-0.92352891) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(-1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975824) q[0];
sx q[0];
rz(-2.0799461) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(-1.7890731) q[1];
sx q[1];
rz(-1.0357772) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3269743) q[0];
sx q[0];
rz(-2.1332392) q[0];
sx q[0];
rz(-0.66280611) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1582583) q[2];
sx q[2];
rz(-0.87450714) q[2];
sx q[2];
rz(-1.0227026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28029728) q[1];
sx q[1];
rz(-1.1251083) q[1];
sx q[1];
rz(2.6764826) q[1];
x q[2];
rz(-0.61061065) q[3];
sx q[3];
rz(-2.7770677) q[3];
sx q[3];
rz(-0.025346905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(2.7374632) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-1.1288246) q[3];
sx q[3];
rz(2.9742441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(-2.6722233) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(-0.00096360047) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199333) q[0];
sx q[0];
rz(-2.1533305) q[0];
sx q[0];
rz(-0.68908738) q[0];
x q[1];
rz(-1.5702815) q[2];
sx q[2];
rz(-1.6948023) q[2];
sx q[2];
rz(2.4027783) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5597792) q[1];
sx q[1];
rz(-2.1372791) q[1];
sx q[1];
rz(0.15484667) q[1];
x q[2];
rz(1.6603052) q[3];
sx q[3];
rz(-2.138391) q[3];
sx q[3];
rz(-1.9247418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.5116189) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(2.3305437) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(-3.1276935) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1902996) q[0];
sx q[0];
rz(-2.2409029) q[0];
sx q[0];
rz(-2.7850889) q[0];
rz(2.6351392) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(-2.8809663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.168419) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(-0.14194685) q[0];
rz(-0.79029681) q[2];
sx q[2];
rz(-0.34966636) q[2];
sx q[2];
rz(-2.5555573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.063559859) q[1];
sx q[1];
rz(-1.5531925) q[1];
sx q[1];
rz(-2.5672008) q[1];
rz(0.77633206) q[3];
sx q[3];
rz(-0.20271248) q[3];
sx q[3];
rz(0.67938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2300718) q[2];
sx q[2];
rz(-1.4804966) q[2];
sx q[2];
rz(-2.9122706) q[2];
rz(-2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.649958) q[0];
rz(-1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(1.414149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0522239) q[0];
sx q[0];
rz(-1.6181237) q[0];
sx q[0];
rz(0.27520673) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51775198) q[2];
sx q[2];
rz(-1.8038097) q[2];
sx q[2];
rz(-1.053996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27281877) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(-2.6150319) q[1];
x q[2];
rz(-1.9023444) q[3];
sx q[3];
rz(-2.9132531) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-0.61855519) q[2];
sx q[2];
rz(0.027739851) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.490482) q[3];
sx q[3];
rz(1.3940575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(2.7391403) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6399022) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(0.34696607) q[0];
x q[1];
rz(-1.5791113) q[2];
sx q[2];
rz(-1.5247716) q[2];
sx q[2];
rz(0.41551057) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9060045) q[1];
sx q[1];
rz(-1.2885805) q[1];
sx q[1];
rz(-1.3868335) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6472858) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(-2.4793712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.014331269) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(-2.4086337) q[0];
rz(-0.14006242) q[1];
sx q[1];
rz(-0.99761325) q[1];
sx q[1];
rz(-1.0345116) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93787545) q[0];
sx q[0];
rz(-0.28521842) q[0];
sx q[0];
rz(-2.2144149) q[0];
rz(-0.10266281) q[2];
sx q[2];
rz(-1.5005939) q[2];
sx q[2];
rz(1.3753124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4701177) q[1];
sx q[1];
rz(-2.5744994) q[1];
sx q[1];
rz(1.5017897) q[1];
rz(-0.23969527) q[3];
sx q[3];
rz(-0.42907676) q[3];
sx q[3];
rz(-0.33340463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2404279) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(2.365716) q[2];
rz(0.72426978) q[3];
sx q[3];
rz(-0.80562076) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-0.7787849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22395615) q[0];
sx q[0];
rz(-0.42376873) q[0];
sx q[0];
rz(0.64026041) q[0];
x q[1];
rz(0.2741371) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(-1.0711087) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89214954) q[1];
sx q[1];
rz(-1.3378694) q[1];
sx q[1];
rz(1.9473142) q[1];
rz(1.7917463) q[3];
sx q[3];
rz(-1.6337992) q[3];
sx q[3];
rz(0.27730478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4497711) q[2];
sx q[2];
rz(-1.3042973) q[2];
sx q[2];
rz(0.71643913) q[2];
rz(-1.4572432) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0338106) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(0.17627136) q[1];
sx q[1];
rz(-1.9468032) q[1];
sx q[1];
rz(2.418628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3176214) q[0];
sx q[0];
rz(-0.184632) q[0];
sx q[0];
rz(-2.1872107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16913551) q[2];
sx q[2];
rz(-1.198487) q[2];
sx q[2];
rz(-1.7793836) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55262676) q[1];
sx q[1];
rz(-2.6347661) q[1];
sx q[1];
rz(-1.456702) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1687752) q[3];
sx q[3];
rz(-0.82377454) q[3];
sx q[3];
rz(-2.0603501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-1.107639) q[3];
sx q[3];
rz(-2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(1.0992959) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(-2.6635086) q[2];
sx q[2];
rz(-2.1838084) q[2];
sx q[2];
rz(2.9391391) q[2];
rz(2.5476534) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
