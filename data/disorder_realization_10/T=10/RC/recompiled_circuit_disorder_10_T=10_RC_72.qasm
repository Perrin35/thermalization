OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878108) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(-1.5009297) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2379797) q[2];
sx q[2];
rz(-1.4822072) q[2];
sx q[2];
rz(-0.36117902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5391985) q[1];
sx q[1];
rz(-1.0326003) q[1];
sx q[1];
rz(-2.5799275) q[1];
x q[2];
rz(-0.88793036) q[3];
sx q[3];
rz(-3.008932) q[3];
sx q[3];
rz(-2.5761029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4937218) q[2];
sx q[2];
rz(-1.9335258) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(-0.99938756) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(2.5527111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988929) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(0.47168628) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6149711) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(-1.4814266) q[0];
rz(0.50498982) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(1.9976975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7029611) q[1];
sx q[1];
rz(-1.8384117) q[1];
sx q[1];
rz(2.5665934) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3470115) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(-0.97298813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.9886118) q[2];
sx q[2];
rz(0.96898752) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(-1.0415174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4330924) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(2.95978) q[0];
rz(-1.1026985) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91711125) q[0];
sx q[0];
rz(-0.79229504) q[0];
sx q[0];
rz(2.2729421) q[0];
rz(-pi) q[1];
rz(-0.45708926) q[2];
sx q[2];
rz(-0.79490137) q[2];
sx q[2];
rz(2.7283816) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1959343) q[1];
sx q[1];
rz(-3.0155026) q[1];
sx q[1];
rz(-3.0581711) q[1];
rz(1.4643747) q[3];
sx q[3];
rz(-0.94821804) q[3];
sx q[3];
rz(-2.4604083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.4870746) q[2];
sx q[2];
rz(-0.96763119) q[2];
rz(0.72757059) q[3];
sx q[3];
rz(-1.8811767) q[3];
sx q[3];
rz(2.9038866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(0.63823429) q[0];
rz(-2.0137285) q[1];
sx q[1];
rz(-0.82740873) q[1];
sx q[1];
rz(1.9086054) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5467984) q[0];
sx q[0];
rz(-1.5017171) q[0];
sx q[0];
rz(-0.49537201) q[0];
rz(-pi) q[1];
rz(-2.4013176) q[2];
sx q[2];
rz(-2.4057655) q[2];
sx q[2];
rz(2.1528113) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14703688) q[1];
sx q[1];
rz(-1.6068659) q[1];
sx q[1];
rz(-1.3922763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87807699) q[3];
sx q[3];
rz(-1.3948166) q[3];
sx q[3];
rz(-0.39719492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91810742) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(0.81400648) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(1.1842747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994613) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(-0.89170757) q[0];
rz(-1.2437598) q[1];
sx q[1];
rz(-1.7638821) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034886995) q[0];
sx q[0];
rz(-1.5505152) q[0];
sx q[0];
rz(3.1253392) q[0];
rz(-pi) q[1];
rz(-1.6264621) q[2];
sx q[2];
rz(-0.60534436) q[2];
sx q[2];
rz(-2.2948613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6032431) q[1];
sx q[1];
rz(-1.0298567) q[1];
sx q[1];
rz(0.69640883) q[1];
x q[2];
rz(1.0605293) q[3];
sx q[3];
rz(-2.3465996) q[3];
sx q[3];
rz(1.0802964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1084958) q[2];
sx q[2];
rz(-2.1704845) q[2];
sx q[2];
rz(-2.6203716) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.228046) q[3];
sx q[3];
rz(-1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8591156) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21340428) q[0];
sx q[0];
rz(-1.8367935) q[0];
sx q[0];
rz(0.022797419) q[0];
rz(-1.5224783) q[2];
sx q[2];
rz(-1.2836576) q[2];
sx q[2];
rz(1.552812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1212335) q[1];
sx q[1];
rz(-1.4660144) q[1];
sx q[1];
rz(-1.1377513) q[1];
rz(-0.029929786) q[3];
sx q[3];
rz(-0.79513351) q[3];
sx q[3];
rz(-3.0628169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1569415) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(-2.8267982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(-0.21952195) q[1];
sx q[1];
rz(-1.6794645) q[1];
sx q[1];
rz(0.25442466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2339904) q[0];
sx q[0];
rz(-1.6798717) q[0];
sx q[0];
rz(1.8926189) q[0];
rz(2.1830325) q[2];
sx q[2];
rz(-1.6974291) q[2];
sx q[2];
rz(2.3938092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88528819) q[1];
sx q[1];
rz(-1.8418152) q[1];
sx q[1];
rz(-2.6726252) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80612225) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(-0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97757942) q[2];
sx q[2];
rz(-1.7501202) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(-1.0036489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106237) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(-2.3176106) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5751858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2306838) q[0];
sx q[0];
rz(-2.3678603) q[0];
sx q[0];
rz(0.45549972) q[0];
x q[1];
rz(0.17922108) q[2];
sx q[2];
rz(-1.0374829) q[2];
sx q[2];
rz(-2.3537677) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7745061) q[1];
sx q[1];
rz(-1.9371867) q[1];
sx q[1];
rz(2.7584502) q[1];
x q[2];
rz(2.4616562) q[3];
sx q[3];
rz(-2.5222062) q[3];
sx q[3];
rz(3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(0.27754647) q[2];
rz(1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(2.1267166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16185109) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(-1.6850527) q[0];
rz(2.5121571) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(2.004752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0658543) q[0];
sx q[0];
rz(-2.542001) q[0];
sx q[0];
rz(1.0435186) q[0];
x q[1];
rz(-2.7669737) q[2];
sx q[2];
rz(-1.1748474) q[2];
sx q[2];
rz(1.7916726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5053619) q[1];
sx q[1];
rz(-2.8103235) q[1];
sx q[1];
rz(-1.2760217) q[1];
rz(0.40739079) q[3];
sx q[3];
rz(-2.3630777) q[3];
sx q[3];
rz(2.9152169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4920766) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(1.0409522) q[2];
rz(3.1395636) q[3];
sx q[3];
rz(-2.7414331) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8390389) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(-2.1955406) q[0];
rz(-0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(0.61202234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6291954) q[0];
sx q[0];
rz(-1.1149659) q[0];
sx q[0];
rz(-0.7645316) q[0];
rz(-pi) q[1];
rz(0.34531784) q[2];
sx q[2];
rz(-1.3441663) q[2];
sx q[2];
rz(-0.40030865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40572383) q[1];
sx q[1];
rz(-1.4475665) q[1];
sx q[1];
rz(2.481639) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53491433) q[3];
sx q[3];
rz(-1.2423008) q[3];
sx q[3];
rz(-1.2319777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(2.4882312) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(-0.70070926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(-2.4621261) q[2];
sx q[2];
rz(-1.1107399) q[2];
sx q[2];
rz(-2.6723292) q[2];
rz(-1.8176953) q[3];
sx q[3];
rz(-0.3331475) q[3];
sx q[3];
rz(-3.0039136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
