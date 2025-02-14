OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.06935057) q[0];
sx q[0];
rz(-2.0877593) q[0];
sx q[0];
rz(-1.2487489) q[0];
rz(-1.2381923) q[1];
sx q[1];
rz(3.5817322) q[1];
sx q[1];
rz(11.323827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71860367) q[0];
sx q[0];
rz(-1.1866633) q[0];
sx q[0];
rz(0.61189135) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79808198) q[2];
sx q[2];
rz(-1.5028364) q[2];
sx q[2];
rz(-0.77347212) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.203071) q[1];
sx q[1];
rz(-1.8932098) q[1];
sx q[1];
rz(1.5713657) q[1];
rz(-0.60939021) q[3];
sx q[3];
rz(-1.8577788) q[3];
sx q[3];
rz(-0.089394102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0480463) q[2];
sx q[2];
rz(-2.8600433) q[2];
sx q[2];
rz(1.2151037) q[2];
rz(1.18527) q[3];
sx q[3];
rz(-0.90648854) q[3];
sx q[3];
rz(1.5989446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95618653) q[0];
sx q[0];
rz(-0.39919272) q[0];
sx q[0];
rz(1.6831552) q[0];
rz(-0.1419119) q[1];
sx q[1];
rz(-1.2071573) q[1];
sx q[1];
rz(-1.9292319) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6120554) q[0];
sx q[0];
rz(-2.3891695) q[0];
sx q[0];
rz(1.0616395) q[0];
x q[1];
rz(-1.5681015) q[2];
sx q[2];
rz(-2.1629253) q[2];
sx q[2];
rz(0.58838974) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72313655) q[1];
sx q[1];
rz(-0.89913705) q[1];
sx q[1];
rz(3.0211012) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47745269) q[3];
sx q[3];
rz(-2.4702854) q[3];
sx q[3];
rz(2.6539913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.045018) q[2];
sx q[2];
rz(-0.18288945) q[2];
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
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69979954) q[0];
sx q[0];
rz(-2.4213591) q[0];
sx q[0];
rz(-2.0752456) q[0];
rz(-2.9329246) q[1];
sx q[1];
rz(-0.51869789) q[1];
sx q[1];
rz(0.64812237) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635669) q[0];
sx q[0];
rz(-1.5180249) q[0];
sx q[0];
rz(0.99356243) q[0];
x q[1];
rz(2.6510973) q[2];
sx q[2];
rz(-1.5083665) q[2];
sx q[2];
rz(1.708781) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0433139) q[1];
sx q[1];
rz(-0.85426353) q[1];
sx q[1];
rz(-1.0745506) q[1];
x q[2];
rz(-2.4345458) q[3];
sx q[3];
rz(-1.505369) q[3];
sx q[3];
rz(-2.9312627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85731792) q[2];
sx q[2];
rz(-1.6782574) q[2];
sx q[2];
rz(0.54764444) q[2];
rz(-2.137843) q[3];
sx q[3];
rz(-0.32310969) q[3];
sx q[3];
rz(-0.16166648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169287) q[0];
sx q[0];
rz(-2.014761) q[0];
sx q[0];
rz(-2.774985) q[0];
rz(-1.2855351) q[1];
sx q[1];
rz(-1.5204241) q[1];
sx q[1];
rz(-2.3191648) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7912589) q[0];
sx q[0];
rz(-0.21203498) q[0];
sx q[0];
rz(-0.37104443) q[0];
rz(-1.0210432) q[2];
sx q[2];
rz(-0.31347358) q[2];
sx q[2];
rz(0.18804729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3318204) q[1];
sx q[1];
rz(-1.0882411) q[1];
sx q[1];
rz(-0.81565522) q[1];
x q[2];
rz(-0.0068914135) q[3];
sx q[3];
rz(-2.8662196) q[3];
sx q[3];
rz(1.7100348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.503868) q[2];
sx q[2];
rz(-0.34054264) q[2];
sx q[2];
rz(1.5857504) q[2];
rz(-1.2268892) q[3];
sx q[3];
rz(-2.5033958) q[3];
sx q[3];
rz(-1.816412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9521088) q[0];
sx q[0];
rz(-1.1186849) q[0];
sx q[0];
rz(-2.191191) q[0];
rz(2.9474126) q[1];
sx q[1];
rz(-1.8870528) q[1];
sx q[1];
rz(-2.8607184) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062894363) q[0];
sx q[0];
rz(-2.4655229) q[0];
sx q[0];
rz(1.990114) q[0];
x q[1];
rz(-1.5834278) q[2];
sx q[2];
rz(-1.4213741) q[2];
sx q[2];
rz(2.0858425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9431587) q[1];
sx q[1];
rz(-1.4869964) q[1];
sx q[1];
rz(-1.6479302) q[1];
rz(-0.30539565) q[3];
sx q[3];
rz(-1.0825233) q[3];
sx q[3];
rz(-2.0528169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6516271) q[2];
sx q[2];
rz(-1.4474892) q[2];
sx q[2];
rz(0.70072407) q[2];
rz(1.0939595) q[3];
sx q[3];
rz(-2.14812) q[3];
sx q[3];
rz(-2.3195364) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1489498) q[0];
sx q[0];
rz(-2.094291) q[0];
sx q[0];
rz(-2.6348422) q[0];
rz(2.0172334) q[1];
sx q[1];
rz(-1.9729112) q[1];
sx q[1];
rz(-2.7728424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4418954) q[0];
sx q[0];
rz(-1.3831733) q[0];
sx q[0];
rz(-0.86098598) q[0];
rz(-pi) q[1];
rz(1.1738335) q[2];
sx q[2];
rz(-0.78437524) q[2];
sx q[2];
rz(-2.337817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3974187) q[1];
sx q[1];
rz(-1.6622006) q[1];
sx q[1];
rz(-1.869702) q[1];
rz(-2.0268129) q[3];
sx q[3];
rz(-1.9590833) q[3];
sx q[3];
rz(-0.20336313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90699482) q[2];
sx q[2];
rz(-0.58595053) q[2];
sx q[2];
rz(-2.9242945) q[2];
rz(-1.6654738) q[3];
sx q[3];
rz(-2.3264591) q[3];
sx q[3];
rz(-2.7998717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23389626) q[0];
sx q[0];
rz(-0.32231575) q[0];
sx q[0];
rz(0.79793683) q[0];
rz(-1.4403053) q[1];
sx q[1];
rz(-2.1013575) q[1];
sx q[1];
rz(-0.61680102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9413404) q[0];
sx q[0];
rz(-1.6464707) q[0];
sx q[0];
rz(-2.5067634) q[0];
x q[1];
rz(-1.2347414) q[2];
sx q[2];
rz(-1.0288887) q[2];
sx q[2];
rz(-1.0934747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7239769) q[1];
sx q[1];
rz(-0.25088746) q[1];
sx q[1];
rz(2.1881359) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2374898) q[3];
sx q[3];
rz(-2.1573503) q[3];
sx q[3];
rz(1.3400638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1282318) q[2];
sx q[2];
rz(-1.1527454) q[2];
sx q[2];
rz(2.4398003) q[2];
rz(-0.93891406) q[3];
sx q[3];
rz(-2.2434668) q[3];
sx q[3];
rz(-1.8963337) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.668648) q[0];
sx q[0];
rz(-2.9845147) q[0];
sx q[0];
rz(0.12338403) q[0];
rz(0.75048796) q[1];
sx q[1];
rz(-1.4742943) q[1];
sx q[1];
rz(-1.0184681) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3833419) q[0];
sx q[0];
rz(-1.4660379) q[0];
sx q[0];
rz(-1.7927756) q[0];
x q[1];
rz(-0.84853402) q[2];
sx q[2];
rz(-0.64788514) q[2];
sx q[2];
rz(2.2415438) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2606517) q[1];
sx q[1];
rz(-0.48704942) q[1];
sx q[1];
rz(-0.080663514) q[1];
rz(-1.4470149) q[3];
sx q[3];
rz(-0.63237337) q[3];
sx q[3];
rz(3.0770709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6259152) q[2];
sx q[2];
rz(-1.6720142) q[2];
sx q[2];
rz(0.53543004) q[2];
rz(-1.911602) q[3];
sx q[3];
rz(-2.1401236) q[3];
sx q[3];
rz(3.1297704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5984421) q[0];
sx q[0];
rz(-2.4813528) q[0];
sx q[0];
rz(1.2793596) q[0];
rz(1.2641501) q[1];
sx q[1];
rz(-1.2707571) q[1];
sx q[1];
rz(2.989891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8063552) q[0];
sx q[0];
rz(-1.6280109) q[0];
sx q[0];
rz(-2.7430659) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13387605) q[2];
sx q[2];
rz(-1.5508442) q[2];
sx q[2];
rz(-0.48568113) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6231767) q[1];
sx q[1];
rz(-2.4881426) q[1];
sx q[1];
rz(1.0423686) q[1];
x q[2];
rz(2.1131198) q[3];
sx q[3];
rz(-0.92307011) q[3];
sx q[3];
rz(-2.5437989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.58574) q[2];
sx q[2];
rz(-2.2057605) q[2];
sx q[2];
rz(1.9110511) q[2];
rz(-1.7963643) q[3];
sx q[3];
rz(-1.3627005) q[3];
sx q[3];
rz(-1.8432553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96974385) q[0];
sx q[0];
rz(-1.3405223) q[0];
sx q[0];
rz(-0.44542435) q[0];
rz(-0.42462665) q[1];
sx q[1];
rz(-1.5698965) q[1];
sx q[1];
rz(0.62579036) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9827756) q[0];
sx q[0];
rz(-2.3390769) q[0];
sx q[0];
rz(-2.3071852) q[0];
rz(1.2568057) q[2];
sx q[2];
rz(-0.50334785) q[2];
sx q[2];
rz(-1.373601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0808472) q[1];
sx q[1];
rz(-2.0909799) q[1];
sx q[1];
rz(2.84723) q[1];
x q[2];
rz(-0.24769737) q[3];
sx q[3];
rz(-0.15532914) q[3];
sx q[3];
rz(-2.4904136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9618535) q[2];
sx q[2];
rz(-1.0280131) q[2];
sx q[2];
rz(1.4614159) q[2];
rz(-3.0840868) q[3];
sx q[3];
rz(-1.4534566) q[3];
sx q[3];
rz(2.1632975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7326603) q[0];
sx q[0];
rz(-1.4807777) q[0];
sx q[0];
rz(-2.5430191) q[0];
rz(2.9231425) q[1];
sx q[1];
rz(-1.5427867) q[1];
sx q[1];
rz(-1.3839518) q[1];
rz(0.54616164) q[2];
sx q[2];
rz(-1.4991789) q[2];
sx q[2];
rz(-3.0747902) q[2];
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
