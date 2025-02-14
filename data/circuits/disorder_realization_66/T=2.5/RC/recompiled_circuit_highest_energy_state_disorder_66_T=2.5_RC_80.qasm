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
rz(2.1386327) q[0];
sx q[0];
rz(-0.47337368) q[0];
sx q[0];
rz(2.0546497) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87693857) q[0];
sx q[0];
rz(-1.042141) q[0];
sx q[0];
rz(2.1318046) q[0];
x q[1];
rz(1.126312) q[2];
sx q[2];
rz(-0.51156564) q[2];
sx q[2];
rz(1.4064521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8352141) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(-1.6914961) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5974724) q[3];
sx q[3];
rz(-1.6427338) q[3];
sx q[3];
rz(-1.6202591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5272556) q[2];
sx q[2];
rz(-1.273512) q[2];
sx q[2];
rz(2.7296076) q[2];
rz(2.8968503) q[3];
sx q[3];
rz(-0.63260806) q[3];
sx q[3];
rz(-2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43993846) q[0];
sx q[0];
rz(-2.169862) q[0];
sx q[0];
rz(-2.7213851) q[0];
rz(0.48928753) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(-1.18321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8156957) q[0];
sx q[0];
rz(-2.1854379) q[0];
sx q[0];
rz(-2.6026898) q[0];
x q[1];
rz(-1.2474044) q[2];
sx q[2];
rz(-1.8696491) q[2];
sx q[2];
rz(-3.1355203) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72231659) q[1];
sx q[1];
rz(-2.4240383) q[1];
sx q[1];
rz(-2.5408387) q[1];
x q[2];
rz(1.39721) q[3];
sx q[3];
rz(-1.6866444) q[3];
sx q[3];
rz(0.66044745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45541397) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(-0.14032826) q[2];
rz(-2.8651107) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8088733) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(1.3746388) q[0];
rz(0.91284347) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(-1.0309781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682273) q[0];
sx q[0];
rz(-1.388167) q[0];
sx q[0];
rz(-1.2398861) q[0];
rz(-pi) q[1];
rz(1.0101843) q[2];
sx q[2];
rz(-1.2282145) q[2];
sx q[2];
rz(0.48784384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.43592855) q[1];
sx q[1];
rz(-1.9829653) q[1];
sx q[1];
rz(0.010910587) q[1];
rz(-pi) q[2];
rz(1.7632906) q[3];
sx q[3];
rz(-1.2519426) q[3];
sx q[3];
rz(-1.1553193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1294407) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(0.63808179) q[2];
rz(0.10073999) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(-2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3070372) q[0];
sx q[0];
rz(-1.5191673) q[0];
sx q[0];
rz(1.7592953) q[0];
rz(2.8221829) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(-0.75739783) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35505193) q[0];
sx q[0];
rz(-3.0664223) q[0];
sx q[0];
rz(3.0274903) q[0];
rz(-0.14230911) q[2];
sx q[2];
rz(-1.6985296) q[2];
sx q[2];
rz(-0.74554986) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74906236) q[1];
sx q[1];
rz(-0.99732256) q[1];
sx q[1];
rz(0.5368781) q[1];
rz(0.60507617) q[3];
sx q[3];
rz(-1.1272484) q[3];
sx q[3];
rz(-3.0435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3697529) q[2];
sx q[2];
rz(-0.87551337) q[2];
sx q[2];
rz(-2.3337951) q[2];
rz(-1.4098903) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(-2.4894449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(-0.98528969) q[0];
rz(-1.4957042) q[1];
sx q[1];
rz(-2.0031877) q[1];
sx q[1];
rz(2.2122502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7791361) q[0];
sx q[0];
rz(-2.1497288) q[0];
sx q[0];
rz(1.8916564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60784908) q[2];
sx q[2];
rz(-1.2298541) q[2];
sx q[2];
rz(0.268595) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73354355) q[1];
sx q[1];
rz(-2.0966623) q[1];
sx q[1];
rz(0.2307363) q[1];
rz(-pi) q[2];
rz(-2.1319904) q[3];
sx q[3];
rz(-1.5586886) q[3];
sx q[3];
rz(2.5502513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.052496584) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(-2.055114) q[2];
rz(0.9440445) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(-2.0199147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97359598) q[0];
sx q[0];
rz(-0.43742988) q[0];
sx q[0];
rz(-0.22015372) q[0];
rz(1.6379697) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(-2.5880609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93388825) q[0];
sx q[0];
rz(-1.3739062) q[0];
sx q[0];
rz(1.6879199) q[0];
x q[1];
rz(-0.36350162) q[2];
sx q[2];
rz(-1.40503) q[2];
sx q[2];
rz(1.7004418) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1268688) q[1];
sx q[1];
rz(-0.8783825) q[1];
sx q[1];
rz(1.6339373) q[1];
rz(-0.27724482) q[3];
sx q[3];
rz(-0.70477761) q[3];
sx q[3];
rz(0.98990209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0452051) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(-2.9333147) q[2];
rz(-2.0221209) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(-2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.646362) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(-1.0736504) q[0];
rz(-0.13058361) q[1];
sx q[1];
rz(-1.5740296) q[1];
sx q[1];
rz(-2.1317587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9511418) q[0];
sx q[0];
rz(-1.1755687) q[0];
sx q[0];
rz(1.7717096) q[0];
x q[1];
rz(-1.4973348) q[2];
sx q[2];
rz(-1.7188004) q[2];
sx q[2];
rz(-0.55390394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4676795) q[1];
sx q[1];
rz(-0.99930489) q[1];
sx q[1];
rz(1.4614146) q[1];
rz(-pi) q[2];
rz(-2.4773682) q[3];
sx q[3];
rz(-1.8598607) q[3];
sx q[3];
rz(2.5632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69663557) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(2.8182287) q[2];
rz(-0.67445406) q[3];
sx q[3];
rz(-2.2489397) q[3];
sx q[3];
rz(0.023580624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0167639) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-2.0523409) q[0];
rz(-0.59843868) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(0.79900297) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625811) q[0];
sx q[0];
rz(-0.41951734) q[0];
sx q[0];
rz(-0.74504344) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4510462) q[2];
sx q[2];
rz(-0.66968067) q[2];
sx q[2];
rz(-2.0872175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5889783) q[1];
sx q[1];
rz(-1.2689021) q[1];
sx q[1];
rz(-1.4149069) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24183065) q[3];
sx q[3];
rz(-1.5655091) q[3];
sx q[3];
rz(1.9092803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66990718) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(1.8890107) q[2];
rz(2.0910828) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893386) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(-0.5156714) q[0];
rz(1.6853261) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(2.8175443) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0607325) q[0];
sx q[0];
rz(-1.3172842) q[0];
sx q[0];
rz(-2.7297165) q[0];
rz(-pi) q[1];
rz(-2.222175) q[2];
sx q[2];
rz(-1.9371261) q[2];
sx q[2];
rz(1.1446143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4676592) q[1];
sx q[1];
rz(-2.0083462) q[1];
sx q[1];
rz(0.15435855) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.079599722) q[3];
sx q[3];
rz(-2.1270424) q[3];
sx q[3];
rz(-1.2885827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9825762) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(0.70875657) q[2];
rz(-2.3580264) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(-1.5172575) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6056972) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(2.9851483) q[0];
rz(-2.5018196) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(-0.99008647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759686) q[0];
sx q[0];
rz(-2.1556615) q[0];
sx q[0];
rz(-2.3192899) q[0];
rz(-pi) q[1];
rz(-1.9205356) q[2];
sx q[2];
rz(-1.4863401) q[2];
sx q[2];
rz(1.3049558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2459979) q[1];
sx q[1];
rz(-2.5709472) q[1];
sx q[1];
rz(-0.24773189) q[1];
rz(-0.57164945) q[3];
sx q[3];
rz(-0.23861966) q[3];
sx q[3];
rz(2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(0.61335316) q[2];
rz(0.26398811) q[3];
sx q[3];
rz(-1.8050906) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(-1.6773979) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(-2.7395372) q[2];
sx q[2];
rz(-1.9378312) q[2];
sx q[2];
rz(2.1294049) q[2];
rz(3.0861978) q[3];
sx q[3];
rz(-2.5221586) q[3];
sx q[3];
rz(0.4172162) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
