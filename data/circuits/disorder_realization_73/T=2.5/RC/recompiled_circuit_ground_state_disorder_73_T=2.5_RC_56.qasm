OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(-2.8889416) q[0];
sx q[0];
rz(1.2055612) q[0];
rz(-1.128101) q[1];
sx q[1];
rz(4.4681273) q[1];
sx q[1];
rz(11.820643) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099961258) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(3.0801386) q[0];
x q[1];
rz(0.22439297) q[2];
sx q[2];
rz(-1.577415) q[2];
sx q[2];
rz(1.4057019) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28994014) q[1];
sx q[1];
rz(-1.2356241) q[1];
sx q[1];
rz(-0.22214684) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0449355) q[3];
sx q[3];
rz(-1.385194) q[3];
sx q[3];
rz(-1.7222278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2241609) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.4707627) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5690174) q[0];
sx q[0];
rz(-1.9439531) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-0.1145656) q[1];
sx q[1];
rz(0.13470185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5935105) q[0];
sx q[0];
rz(-2.5786886) q[0];
sx q[0];
rz(2.4538732) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7472166) q[2];
sx q[2];
rz(-0.069709502) q[2];
sx q[2];
rz(1.4331872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33011064) q[1];
sx q[1];
rz(-2.5977784) q[1];
sx q[1];
rz(-2.5303929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7487583) q[3];
sx q[3];
rz(-1.4918394) q[3];
sx q[3];
rz(2.3633901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(-2.9888195) q[2];
rz(-1.7759391) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(0.16807817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.633054) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(-2.659944) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.7748723) q[1];
sx q[1];
rz(2.1462323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5115968) q[0];
sx q[0];
rz(-0.42034402) q[0];
sx q[0];
rz(0.98031251) q[0];
rz(-pi) q[1];
rz(1.5103417) q[2];
sx q[2];
rz(-1.619307) q[2];
sx q[2];
rz(-1.2918351) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9650594) q[1];
sx q[1];
rz(-0.73672026) q[1];
sx q[1];
rz(2.6882369) q[1];
x q[2];
rz(-2.1667492) q[3];
sx q[3];
rz(-1.5961093) q[3];
sx q[3];
rz(1.635575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2153726) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(3.0769297) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30508405) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-2.3205561) q[0];
rz(3.0631284) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-0.99748126) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0745463) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(1.9590098) q[0];
x q[1];
rz(1.5115963) q[2];
sx q[2];
rz(-1.534605) q[2];
sx q[2];
rz(2.9179887) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23826829) q[1];
sx q[1];
rz(-1.0449755) q[1];
sx q[1];
rz(-2.8366798) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6281582) q[3];
sx q[3];
rz(-1.093075) q[3];
sx q[3];
rz(-2.3582332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(-1.1537665) q[2];
rz(0.18618259) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(0.33128273) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1368197) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.9983043) q[0];
rz(2.000287) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(-2.6236261) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9310936) q[0];
sx q[0];
rz(-0.8862952) q[0];
sx q[0];
rz(2.4263591) q[0];
rz(-1.5094366) q[2];
sx q[2];
rz(-0.01505919) q[2];
sx q[2];
rz(-1.7565774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3822381) q[1];
sx q[1];
rz(-1.2975946) q[1];
sx q[1];
rz(-1.9135936) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.84359) q[3];
sx q[3];
rz(-1.4299222) q[3];
sx q[3];
rz(-2.6948787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(-0.32773584) q[2];
rz(1.94708) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(0.56185454) q[0];
rz(1.4909164) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(3.0432826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54360169) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(1.4844358) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5709236) q[2];
sx q[2];
rz(-1.5701558) q[2];
sx q[2];
rz(2.8252779) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5206388) q[1];
sx q[1];
rz(-0.47157447) q[1];
sx q[1];
rz(-3.0388799) q[1];
rz(2.6410511) q[3];
sx q[3];
rz(-1.1574452) q[3];
sx q[3];
rz(-2.7374637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(1.0657715) q[2];
rz(2.7417475) q[3];
sx q[3];
rz(-0.53374922) q[3];
sx q[3];
rz(-1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-0.081900224) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(2.0211925) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-0.33946005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0507815) q[0];
sx q[0];
rz(-0.078074038) q[0];
sx q[0];
rz(2.9705621) q[0];
rz(-pi) q[1];
rz(-0.022384251) q[2];
sx q[2];
rz(-2.5867992) q[2];
sx q[2];
rz(0.024352976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.124954) q[1];
sx q[1];
rz(-1.5724239) q[1];
sx q[1];
rz(-1.485199) q[1];
rz(-pi) q[2];
rz(2.8027995) q[3];
sx q[3];
rz(-1.6648035) q[3];
sx q[3];
rz(0.28123048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7772943) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61984396) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(0.92292619) q[0];
rz(-1.6587616) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-3.1373851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18989604) q[0];
sx q[0];
rz(-1.0605264) q[0];
sx q[0];
rz(-2.4416591) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9816859) q[2];
sx q[2];
rz(-1.5762323) q[2];
sx q[2];
rz(1.5856575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1394704) q[1];
sx q[1];
rz(-1.0957901) q[1];
sx q[1];
rz(0.072474555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5291653) q[3];
sx q[3];
rz(-1.53293) q[3];
sx q[3];
rz(1.7671536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5795472) q[2];
sx q[2];
rz(-1.6041218) q[2];
sx q[2];
rz(-1.2008249) q[2];
rz(1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(-2.4414731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-2.6926079) q[0];
sx q[0];
rz(-1.2196983) q[0];
rz(1.3297184) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(-0.16389287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71918286) q[0];
sx q[0];
rz(-1.9253732) q[0];
sx q[0];
rz(2.9585881) q[0];
x q[1];
rz(1.8437949) q[2];
sx q[2];
rz(-2.1381209) q[2];
sx q[2];
rz(1.3162055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6299705) q[1];
sx q[1];
rz(-3.0093806) q[1];
sx q[1];
rz(-1.1019143) q[1];
rz(-pi) q[2];
rz(0.90198646) q[3];
sx q[3];
rz(-0.94873442) q[3];
sx q[3];
rz(-0.40611503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-3.0516477) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(-3.0981787) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(-2.6553335) q[0];
rz(-2.4468415) q[1];
sx q[1];
rz(-2.7943352) q[1];
sx q[1];
rz(-2.6659226) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566475) q[0];
sx q[0];
rz(-1.534919) q[0];
sx q[0];
rz(-3.1079124) q[0];
rz(1.6108405) q[2];
sx q[2];
rz(-0.70654987) q[2];
sx q[2];
rz(-1.5020811) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039363843) q[1];
sx q[1];
rz(-2.1783531) q[1];
sx q[1];
rz(-1.5786912) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6127897) q[3];
sx q[3];
rz(-2.7861054) q[3];
sx q[3];
rz(-2.5606009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67705578) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(1.2976868) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-2.5864351) q[3];
sx q[3];
rz(-0.36778522) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0155335) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(2.2592648) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(-2.7693314) q[2];
sx q[2];
rz(-3.0426171) q[2];
sx q[2];
rz(-0.097886861) q[2];
rz(-1.5658548) q[3];
sx q[3];
rz(-1.3316122) q[3];
sx q[3];
rz(0.022461654) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
