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
rz(-1.128101) q[1];
sx q[1];
rz(-1.815058) q[1];
sx q[1];
rz(2.3958652) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6572031) q[0];
sx q[0];
rz(-1.6307388) q[0];
sx q[0];
rz(-1.7933229) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9171997) q[2];
sx q[2];
rz(-1.5641777) q[2];
sx q[2];
rz(1.7358908) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8516525) q[1];
sx q[1];
rz(-1.2356241) q[1];
sx q[1];
rz(-2.9194458) q[1];
rz(1.9609591) q[3];
sx q[3];
rz(-0.50658617) q[3];
sx q[3];
rz(-2.6449639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91743177) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(-1.4707627) q[2];
rz(1.6338232) q[3];
sx q[3];
rz(-0.016851146) q[3];
sx q[3];
rz(-2.2182218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(1.5731328) q[0];
rz(-0.16954999) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(-3.0068908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5935105) q[0];
sx q[0];
rz(-0.56290409) q[0];
sx q[0];
rz(-2.4538732) q[0];
x q[1];
rz(1.3943761) q[2];
sx q[2];
rz(-3.0718832) q[2];
sx q[2];
rz(-1.4331872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0162279) q[1];
sx q[1];
rz(-2.008359) q[1];
sx q[1];
rz(1.9047649) q[1];
rz(1.7487583) q[3];
sx q[3];
rz(-1.6497532) q[3];
sx q[3];
rz(-2.3633901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(2.9888195) q[2];
rz(1.7759391) q[3];
sx q[3];
rz(-0.036866166) q[3];
sx q[3];
rz(2.9735145) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.633054) q[0];
sx q[0];
rz(-2.3336053) q[0];
sx q[0];
rz(2.659944) q[0];
rz(-0.18474361) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(-0.99536037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1448875) q[0];
sx q[0];
rz(-1.9166244) q[0];
sx q[0];
rz(0.24390999) q[0];
rz(-1.6312509) q[2];
sx q[2];
rz(-1.619307) q[2];
sx q[2];
rz(1.8497576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9650594) q[1];
sx q[1];
rz(-0.73672026) q[1];
sx q[1];
rz(2.6882369) q[1];
x q[2];
rz(-1.5257201) q[3];
sx q[3];
rz(-0.59642506) q[3];
sx q[3];
rz(0.10208043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92622009) q[2];
sx q[2];
rz(-3.0871349) q[2];
sx q[2];
rz(3.0769297) q[2];
rz(-1.0489382) q[3];
sx q[3];
rz(-3.1147396) q[3];
sx q[3];
rz(-1.2803199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8365086) q[0];
sx q[0];
rz(-0.11798141) q[0];
sx q[0];
rz(2.3205561) q[0];
rz(0.07846421) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-2.1441114) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0745463) q[0];
sx q[0];
rz(-0.95854488) q[0];
sx q[0];
rz(-1.1825829) q[0];
rz(-1.0216265) q[2];
sx q[2];
rz(-3.0722174) q[2];
sx q[2];
rz(0.79909426) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8197937) q[1];
sx q[1];
rz(-0.60056486) q[1];
sx q[1];
rz(1.0933881) q[1];
rz(-pi) q[2];
rz(1.6281582) q[3];
sx q[3];
rz(-2.0485176) q[3];
sx q[3];
rz(-0.78335947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0492101) q[2];
sx q[2];
rz(-3.1123078) q[2];
sx q[2];
rz(1.1537665) q[2];
rz(-0.18618259) q[3];
sx q[3];
rz(-0.081705339) q[3];
sx q[3];
rz(0.33128273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004772923) q[0];
sx q[0];
rz(-2.3240219) q[0];
sx q[0];
rz(1.1432884) q[0];
rz(2.000287) q[1];
sx q[1];
rz(-0.80991304) q[1];
sx q[1];
rz(-0.51796651) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2835613) q[0];
sx q[0];
rz(-2.1037344) q[0];
sx q[0];
rz(-2.3951247) q[0];
rz(-pi) q[1];
rz(-0.00092351726) q[2];
sx q[2];
rz(-1.5557655) q[2];
sx q[2];
rz(1.3236486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.458788) q[1];
sx q[1];
rz(-0.43495754) q[1];
sx q[1];
rz(-0.87587237) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4235184) q[3];
sx q[3];
rz(-1.8657576) q[3];
sx q[3];
rz(1.9744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3190069) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(2.8138568) q[2];
rz(-1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(0.85668844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(-2.5797381) q[0];
rz(-1.4909164) q[1];
sx q[1];
rz(-1.5166538) q[1];
sx q[1];
rz(3.0432826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.597991) q[0];
sx q[0];
rz(-2.7181135) q[0];
sx q[0];
rz(-1.6571568) q[0];
rz(2.9454284) q[2];
sx q[2];
rz(-3.1409396) q[2];
sx q[2];
rz(-0.12015039) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1001818) q[1];
sx q[1];
rz(-1.5242002) q[1];
sx q[1];
rz(-0.46943922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0344072) q[3];
sx q[3];
rz(-1.1157728) q[3];
sx q[3];
rz(-2.1912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.085122846) q[2];
sx q[2];
rz(-0.19530185) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(-1.8736418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.6837233) q[0];
rz(-1.1204002) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-0.33946005) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65050478) q[0];
sx q[0];
rz(-1.5840713) q[0];
sx q[0];
rz(3.0646532) q[0];
x q[1];
rz(-1.5846662) q[2];
sx q[2];
rz(-1.0161581) q[2];
sx q[2];
rz(-0.0019794606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.554018) q[1];
sx q[1];
rz(-1.4851991) q[1];
sx q[1];
rz(3.1399591) q[1];
x q[2];
rz(-0.27642997) q[3];
sx q[3];
rz(-2.7904841) q[3];
sx q[3];
rz(2.1123667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3642984) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(-0.36336362) q[2];
rz(-1.0904788) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(2.0232078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.61984396) q[0];
sx q[0];
rz(-2.813297) q[0];
sx q[0];
rz(2.2186665) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(0.0042075687) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370458) q[0];
sx q[0];
rz(-0.97386375) q[0];
sx q[0];
rz(2.202522) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1356631) q[2];
sx q[2];
rz(-1.9816795) q[2];
sx q[2];
rz(0.012492736) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.15531047) q[1];
sx q[1];
rz(-2.661507) q[1];
sx q[1];
rz(-1.4309149) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0757853) q[3];
sx q[3];
rz(-0.61344693) q[3];
sx q[3];
rz(-2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5795472) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.2008249) q[2];
rz(-1.7480525) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30018184) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(-1.9218943) q[0];
rz(-1.3297184) q[1];
sx q[1];
rz(-1.1390353) q[1];
sx q[1];
rz(0.16389287) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4224098) q[0];
sx q[0];
rz(-1.9253732) q[0];
sx q[0];
rz(2.9585881) q[0];
rz(-pi) q[1];
rz(-1.8437949) q[2];
sx q[2];
rz(-1.0034717) q[2];
sx q[2];
rz(-1.8253872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6170609) q[1];
sx q[1];
rz(-1.5111898) q[1];
sx q[1];
rz(-1.6888794) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4285942) q[3];
sx q[3];
rz(-0.87942356) q[3];
sx q[3];
rz(0.52934968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71358744) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(-0.62291992) q[2];
rz(0.04341393) q[3];
sx q[3];
rz(-0.89689887) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6503705) q[0];
sx q[0];
rz(-0.0064042052) q[0];
sx q[0];
rz(0.48625913) q[0];
rz(0.69475118) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(-0.4756701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18494517) q[0];
sx q[0];
rz(-1.534919) q[0];
sx q[0];
rz(-0.033680276) q[0];
rz(-1.6108405) q[2];
sx q[2];
rz(-2.4350428) q[2];
sx q[2];
rz(-1.5020811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5359395) q[1];
sx q[1];
rz(-1.5772784) q[1];
sx q[1];
rz(2.5340213) q[1];
x q[2];
rz(2.6127897) q[3];
sx q[3];
rz(-2.7861054) q[3];
sx q[3];
rz(0.58099174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4645369) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(2.2592648) q[1];
sx q[1];
rz(-0.066451646) q[1];
sx q[1];
rz(0.8393504) q[1];
rz(-0.37226128) q[2];
sx q[2];
rz(-0.098975565) q[2];
sx q[2];
rz(3.0437058) q[2];
rz(2.9024057) q[3];
sx q[3];
rz(-1.5659955) q[3];
sx q[3];
rz(-1.5495054) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
