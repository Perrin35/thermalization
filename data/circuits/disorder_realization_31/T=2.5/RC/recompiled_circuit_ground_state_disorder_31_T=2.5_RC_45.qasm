OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0722421) q[0];
sx q[0];
rz(-1.0538333) q[0];
sx q[0];
rz(-1.8928438) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(-2.7014531) q[1];
sx q[1];
rz(1.8990489) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71860367) q[0];
sx q[0];
rz(-1.1866633) q[0];
sx q[0];
rz(2.5297013) q[0];
rz(-pi) q[1];
rz(0.79808198) q[2];
sx q[2];
rz(-1.6387562) q[2];
sx q[2];
rz(0.77347212) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63245508) q[1];
sx q[1];
rz(-1.5702562) q[1];
sx q[1];
rz(-2.8191791) q[1];
rz(-pi) q[2];
rz(1.2253157) q[3];
sx q[3];
rz(-0.98967797) q[3];
sx q[3];
rz(1.2863359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(1.9264889) q[2];
rz(1.18527) q[3];
sx q[3];
rz(-2.2351041) q[3];
sx q[3];
rz(-1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95618653) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(-1.6831552) q[0];
rz(0.1419119) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(1.9292319) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0182918) q[0];
sx q[0];
rz(-2.2102076) q[0];
sx q[0];
rz(-2.713504) q[0];
rz(-pi) q[1];
rz(-1.5681015) q[2];
sx q[2];
rz(-0.97866733) q[2];
sx q[2];
rz(2.5532029) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77245678) q[1];
sx q[1];
rz(-1.6650272) q[1];
sx q[1];
rz(2.2460031) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5271537) q[3];
sx q[3];
rz(-1.2809291) q[3];
sx q[3];
rz(1.4680999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.045018) q[2];
sx q[2];
rz(-2.9587032) q[2];
sx q[2];
rz(2.1796687) q[2];
rz(2.9436881) q[3];
sx q[3];
rz(-1.3062545) q[3];
sx q[3];
rz(-1.643868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69979954) q[0];
sx q[0];
rz(-0.72023359) q[0];
sx q[0];
rz(1.066347) q[0];
rz(-0.20866808) q[1];
sx q[1];
rz(-0.51869789) q[1];
sx q[1];
rz(2.4934703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3828165) q[0];
sx q[0];
rz(-0.99446873) q[0];
sx q[0];
rz(3.0786425) q[0];
rz(0.13193138) q[2];
sx q[2];
rz(-0.49413097) q[2];
sx q[2];
rz(3.1199093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40884205) q[1];
sx q[1];
rz(-2.2957845) q[1];
sx q[1];
rz(-0.50028657) q[1];
rz(0.1005248) q[3];
sx q[3];
rz(-2.4320452) q[3];
sx q[3];
rz(-1.8574992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85731792) q[2];
sx q[2];
rz(-1.4633353) q[2];
sx q[2];
rz(0.54764444) q[2];
rz(2.137843) q[3];
sx q[3];
rz(-0.32310969) q[3];
sx q[3];
rz(-2.9799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52466398) q[0];
sx q[0];
rz(-2.014761) q[0];
sx q[0];
rz(2.774985) q[0];
rz(1.2855351) q[1];
sx q[1];
rz(-1.6211685) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4124968) q[0];
sx q[0];
rz(-1.7682045) q[0];
sx q[0];
rz(-1.4928994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16777487) q[2];
sx q[2];
rz(-1.3047403) q[2];
sx q[2];
rz(0.38420907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3144223) q[1];
sx q[1];
rz(-2.2232375) q[1];
sx q[1];
rz(0.62364044) q[1];
rz(1.5688492) q[3];
sx q[3];
rz(-1.8461627) q[3];
sx q[3];
rz(1.438719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.503868) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(-1.5857504) q[2];
rz(-1.9147035) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(-1.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1894839) q[0];
sx q[0];
rz(-2.0229078) q[0];
sx q[0];
rz(-0.9504016) q[0];
rz(2.9474126) q[1];
sx q[1];
rz(-1.2545398) q[1];
sx q[1];
rz(-0.28087428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9683198) q[0];
sx q[0];
rz(-1.313198) q[0];
sx q[0];
rz(0.93846847) q[0];
rz(-pi) q[1];
rz(-2.9921586) q[2];
sx q[2];
rz(-1.5583056) q[2];
sx q[2];
rz(0.51692671) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54714006) q[1];
sx q[1];
rz(-3.0277589) q[1];
sx q[1];
rz(2.3992541) q[1];
x q[2];
rz(2.0858795) q[3];
sx q[3];
rz(-0.56929811) q[3];
sx q[3];
rz(1.6804939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6516271) q[2];
sx q[2];
rz(-1.6941035) q[2];
sx q[2];
rz(-2.4408686) q[2];
rz(-2.0476332) q[3];
sx q[3];
rz(-0.9934727) q[3];
sx q[3];
rz(-0.82205621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99264282) q[0];
sx q[0];
rz(-1.0473017) q[0];
sx q[0];
rz(-2.6348422) q[0];
rz(2.0172334) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(-2.7728424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0567199) q[0];
sx q[0];
rz(-2.4115926) q[0];
sx q[0];
rz(1.2873136) q[0];
x q[1];
rz(2.314744) q[2];
sx q[2];
rz(-1.294181) q[2];
sx q[2];
rz(-2.6630304) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79850393) q[1];
sx q[1];
rz(-1.868416) q[1];
sx q[1];
rz(0.095620015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1147797) q[3];
sx q[3];
rz(-1.9590833) q[3];
sx q[3];
rz(-0.20336313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90699482) q[2];
sx q[2];
rz(-2.5556421) q[2];
sx q[2];
rz(-0.21729812) q[2];
rz(-1.6654738) q[3];
sx q[3];
rz(-0.81513351) q[3];
sx q[3];
rz(-0.34172094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9076964) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(0.79793683) q[0];
rz(1.7012874) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(-0.61680102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154222) q[0];
sx q[0];
rz(-2.2035193) q[0];
sx q[0];
rz(1.4769082) q[0];
rz(2.6404713) q[2];
sx q[2];
rz(-0.62868147) q[2];
sx q[2];
rz(1.4528265) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.7239769) q[1];
sx q[1];
rz(-2.8907052) q[1];
sx q[1];
rz(2.1881359) q[1];
rz(-pi) q[2];
rz(2.439626) q[3];
sx q[3];
rz(-2.1118374) q[3];
sx q[3];
rz(0.64149414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1282318) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(2.4398003) q[2];
rz(0.93891406) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(-1.2452589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4729446) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(3.0182086) q[0];
rz(-2.3911047) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(-1.0184681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3833419) q[0];
sx q[0];
rz(-1.6755548) q[0];
sx q[0];
rz(1.3488171) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6776601) q[2];
sx q[2];
rz(-2.0407157) q[2];
sx q[2];
rz(1.4063094) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5230549) q[1];
sx q[1];
rz(-1.5330761) q[1];
sx q[1];
rz(-0.48570363) q[1];
rz(-pi) q[2];
rz(-1.4470149) q[3];
sx q[3];
rz(-2.5092193) q[3];
sx q[3];
rz(0.06452175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5156775) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(2.6061626) q[2];
rz(1.911602) q[3];
sx q[3];
rz(-2.1401236) q[3];
sx q[3];
rz(-3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5984421) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(-1.2793596) q[0];
rz(1.2641501) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(2.989891) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3352375) q[0];
sx q[0];
rz(-1.5135817) q[0];
sx q[0];
rz(-2.7430659) q[0];
x q[1];
rz(-0.13387605) q[2];
sx q[2];
rz(-1.5508442) q[2];
sx q[2];
rz(2.6559115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.38167414) q[1];
sx q[1];
rz(-1.8823138) q[1];
sx q[1];
rz(-2.1550234) q[1];
rz(-pi) q[2];
rz(-0.72355481) q[3];
sx q[3];
rz(-1.9949942) q[3];
sx q[3];
rz(-1.8198609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.58574) q[2];
sx q[2];
rz(-2.2057605) q[2];
sx q[2];
rz(-1.2305416) q[2];
rz(-1.3452283) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(1.8432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96974385) q[0];
sx q[0];
rz(-1.8010704) q[0];
sx q[0];
rz(2.6961683) q[0];
rz(-0.42462665) q[1];
sx q[1];
rz(-1.5716962) q[1];
sx q[1];
rz(2.5158023) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75794894) q[0];
sx q[0];
rz(-1.008908) q[0];
sx q[0];
rz(-2.5342219) q[0];
x q[1];
rz(2.9731304) q[2];
sx q[2];
rz(-1.0941774) q[2];
sx q[2];
rz(-2.1229975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6546331) q[1];
sx q[1];
rz(-0.59096293) q[1];
sx q[1];
rz(2.0396359) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8938953) q[3];
sx q[3];
rz(-2.9862635) q[3];
sx q[3];
rz(2.4904136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-2.1135795) q[2];
sx q[2];
rz(-1.4614159) q[2];
rz(-0.05750582) q[3];
sx q[3];
rz(-1.688136) q[3];
sx q[3];
rz(2.1632975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4089324) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(2.9231425) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(-0.54616164) q[2];
sx q[2];
rz(-1.6424137) q[2];
sx q[2];
rz(0.066802468) q[2];
rz(1.5254088) q[3];
sx q[3];
rz(-0.78730351) q[3];
sx q[3];
rz(-1.3135943) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
