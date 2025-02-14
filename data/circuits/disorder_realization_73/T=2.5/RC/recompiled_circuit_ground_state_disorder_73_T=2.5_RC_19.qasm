OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8483228) q[0];
sx q[0];
rz(-0.25265101) q[0];
sx q[0];
rz(-1.2055612) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(-2.3958652) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17240758) q[0];
sx q[0];
rz(-0.23032941) q[0];
sx q[0];
rz(-1.8363097) q[0];
rz(-pi) q[1];
rz(0.22439297) q[2];
sx q[2];
rz(-1.5641777) q[2];
sx q[2];
rz(1.7358908) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9348976) q[1];
sx q[1];
rz(-1.7803915) q[1];
sx q[1];
rz(1.9137726) q[1];
x q[2];
rz(-1.9609591) q[3];
sx q[3];
rz(-2.6350065) q[3];
sx q[3];
rz(0.49662874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.4707627) q[2];
rz(-1.5077695) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725752) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(3.0068908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54808211) q[0];
sx q[0];
rz(-0.56290409) q[0];
sx q[0];
rz(2.4538732) q[0];
rz(1.3943761) q[2];
sx q[2];
rz(-0.069709502) q[2];
sx q[2];
rz(1.4331872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.811482) q[1];
sx q[1];
rz(-0.54381424) q[1];
sx q[1];
rz(-0.61119975) q[1];
x q[2];
rz(1.7487583) q[3];
sx q[3];
rz(-1.4918394) q[3];
sx q[3];
rz(2.3633901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1255101) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(-2.9888195) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(-0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50853866) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(2.659944) q[0];
rz(-0.18474361) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(0.99536037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6516614) q[0];
sx q[0];
rz(-1.341594) q[0];
sx q[0];
rz(-1.9263173) q[0];
rz(-pi) q[1];
rz(0.048599343) q[2];
sx q[2];
rz(-1.5104129) q[2];
sx q[2];
rz(-0.27602613) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0936443) q[1];
sx q[1];
rz(-1.8694832) q[1];
sx q[1];
rz(0.68409749) q[1];
rz(-pi) q[2];
rz(-0.97484346) q[3];
sx q[3];
rz(-1.5454834) q[3];
sx q[3];
rz(-1.5060177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92622009) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(0.064662956) q[2];
rz(-1.0489382) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(-1.8612727) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30508405) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(2.3205561) q[0];
rz(-0.07846421) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-0.99748126) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8686912) q[0];
sx q[0];
rz(-1.2558381) q[0];
sx q[0];
rz(0.64906831) q[0];
rz(3.1053379) q[2];
sx q[2];
rz(-1.5116351) q[2];
sx q[2];
rz(-1.349337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9033244) q[1];
sx q[1];
rz(-2.0966171) q[1];
sx q[1];
rz(-0.30491288) q[1];
rz(-1.5134345) q[3];
sx q[3];
rz(-2.0485176) q[3];
sx q[3];
rz(2.3582332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(-1.9878261) q[2];
rz(-0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(-0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004772923) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.1432884) q[0];
rz(1.1413057) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(2.6236261) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1521027) q[0];
sx q[0];
rz(-0.94606646) q[0];
sx q[0];
rz(-0.89390198) q[0];
rz(-pi) q[1];
rz(-1.5094366) q[2];
sx q[2];
rz(-3.1265335) q[2];
sx q[2];
rz(-1.3850152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3822381) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(-1.9135936) q[1];
rz(-pi) q[2];
rz(1.4235184) q[3];
sx q[3];
rz(-1.8657576) q[3];
sx q[3];
rz(-1.9744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8225857) q[2];
sx q[2];
rz(-2.1876882) q[2];
sx q[2];
rz(-0.32773584) q[2];
rz(1.1945126) q[3];
sx q[3];
rz(-3.0142398) q[3];
sx q[3];
rz(0.85668844) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-2.5797381) q[0];
rz(1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(-3.0432826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033103) q[0];
sx q[0];
rz(-1.9925963) q[0];
sx q[0];
rz(-0.038859239) q[0];
rz(-1.570669) q[2];
sx q[2];
rz(-1.5714368) q[2];
sx q[2];
rz(-0.31631472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.505762) q[1];
sx q[1];
rz(-1.1019076) q[1];
sx q[1];
rz(1.5185578) q[1];
rz(-pi) q[2];
rz(-0.50054153) q[3];
sx q[3];
rz(-1.1574452) q[3];
sx q[3];
rz(0.40412892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(-1.1204002) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-0.33946005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65050478) q[0];
sx q[0];
rz(-1.5840713) q[0];
sx q[0];
rz(-3.0646532) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5869114) q[2];
sx q[2];
rz(-1.5825869) q[2];
sx q[2];
rz(-1.5654711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.124954) q[1];
sx q[1];
rz(-1.5691688) q[1];
sx q[1];
rz(-1.485199) q[1];
rz(1.4711597) q[3];
sx q[3];
rz(-1.2335586) q[3];
sx q[3];
rz(-1.3226313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3642984) q[2];
sx q[2];
rz(-0.0019145049) q[2];
sx q[2];
rz(-0.36336362) q[2];
rz(2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61984396) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(-0.92292619) q[0];
rz(-1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-0.0042075687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85535964) q[0];
sx q[0];
rz(-2.3015733) q[0];
sx q[0];
rz(-0.71536163) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1599067) q[2];
sx q[2];
rz(-1.5653603) q[2];
sx q[2];
rz(1.5559352) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5397268) q[1];
sx q[1];
rz(-1.5063573) q[1];
sx q[1];
rz(-2.0468725) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0757853) q[3];
sx q[3];
rz(-2.5281457) q[3];
sx q[3];
rz(-2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5795472) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(1.9218943) q[0];
rz(-1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.208928) q[0];
sx q[0];
rz(-0.39723662) q[0];
sx q[0];
rz(2.0276638) q[0];
rz(-pi) q[1];
rz(1.8437949) q[2];
sx q[2];
rz(-1.0034717) q[2];
sx q[2];
rz(-1.3162055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6299705) q[1];
sx q[1];
rz(-3.0093806) q[1];
sx q[1];
rz(1.1019143) q[1];
rz(-pi) q[2];
rz(-2.2396062) q[3];
sx q[3];
rz(-0.94873442) q[3];
sx q[3];
rz(-0.40611503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-2.5186727) q[2];
rz(-3.0981787) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(0.77891946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(-2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(-2.6659226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93906389) q[0];
sx q[0];
rz(-3.0923884) q[0];
sx q[0];
rz(-0.81728191) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5307521) q[2];
sx q[2];
rz(-0.70654987) q[2];
sx q[2];
rz(-1.6395115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.025534677) q[1];
sx q[1];
rz(-0.60760159) q[1];
sx q[1];
rz(3.1302384) q[1];
x q[2];
rz(0.52880295) q[3];
sx q[3];
rz(-0.35548726) q[3];
sx q[3];
rz(-2.5606009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67705578) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(-1.2976868) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(-0.36778522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-0.88232782) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(0.092236237) q[2];
sx q[2];
rz(-1.606745) q[2];
sx q[2];
rz(1.8435115) q[2];
rz(-3.121331) q[3];
sx q[3];
rz(-2.9023585) q[3];
sx q[3];
rz(0.0016061847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
