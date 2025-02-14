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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000769) q[0];
sx q[0];
rz(-1.0934663) q[0];
sx q[0];
rz(2.5377089) q[0];
rz(-pi) q[1];
rz(-0.23687266) q[2];
sx q[2];
rz(-2.028596) q[2];
sx q[2];
rz(-1.2352236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30637851) q[1];
sx q[1];
rz(-1.467538) q[1];
sx q[1];
rz(-1.6914961) q[1];
rz(-1.5441203) q[3];
sx q[3];
rz(-1.6427338) q[3];
sx q[3];
rz(1.5213336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5272556) q[2];
sx q[2];
rz(-1.273512) q[2];
sx q[2];
rz(0.41198507) q[2];
rz(2.8968503) q[3];
sx q[3];
rz(-0.63260806) q[3];
sx q[3];
rz(0.28732029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43993846) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(0.42020759) q[0];
rz(2.6523051) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(-1.9583826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32589697) q[0];
sx q[0];
rz(-2.1854379) q[0];
sx q[0];
rz(-0.53890284) q[0];
rz(2.3406896) q[2];
sx q[2];
rz(-2.704853) q[2];
sx q[2];
rz(-2.2855121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1258613) q[1];
sx q[1];
rz(-2.1441048) q[1];
sx q[1];
rz(1.1124951) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11759956) q[3];
sx q[3];
rz(-1.7432074) q[3];
sx q[3];
rz(-2.2109779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6861787) q[2];
sx q[2];
rz(-1.6697465) q[2];
sx q[2];
rz(0.14032826) q[2];
rz(2.8651107) q[3];
sx q[3];
rz(-2.3639207) q[3];
sx q[3];
rz(-0.28237835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33271933) q[0];
sx q[0];
rz(-2.7820899) q[0];
sx q[0];
rz(-1.7669539) q[0];
rz(0.91284347) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(-1.0309781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08374005) q[0];
sx q[0];
rz(-2.7652604) q[0];
sx q[0];
rz(1.0539088) q[0];
rz(-2.7430277) q[2];
sx q[2];
rz(-1.0462648) q[2];
sx q[2];
rz(2.2664859) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46315868) q[1];
sx q[1];
rz(-2.7292876) q[1];
sx q[1];
rz(1.5458471) q[1];
rz(0.52521962) q[3];
sx q[3];
rz(-0.37074836) q[3];
sx q[3];
rz(-1.4300089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-0.63808179) q[2];
rz(3.0408527) q[3];
sx q[3];
rz(-2.1295857) q[3];
sx q[3];
rz(-2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070372) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(1.7592953) q[0];
rz(-2.8221829) q[1];
sx q[1];
rz(-1.6326135) q[1];
sx q[1];
rz(2.3841948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35505193) q[0];
sx q[0];
rz(-0.075170366) q[0];
sx q[0];
rz(0.11410232) q[0];
rz(0.14230911) q[2];
sx q[2];
rz(-1.443063) q[2];
sx q[2];
rz(2.3960428) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74906236) q[1];
sx q[1];
rz(-2.1442701) q[1];
sx q[1];
rz(2.6047146) q[1];
x q[2];
rz(0.60507617) q[3];
sx q[3];
rz(-1.1272484) q[3];
sx q[3];
rz(-3.0435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77183977) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(-2.3337951) q[2];
rz(1.7317023) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0659502) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(2.156303) q[0];
rz(-1.4957042) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(-2.2122502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083223) q[0];
sx q[0];
rz(-2.4886971) q[0];
sx q[0];
rz(0.44954957) q[0];
x q[1];
rz(-2.5856951) q[2];
sx q[2];
rz(-0.68624845) q[2];
sx q[2];
rz(-0.85418073) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73354355) q[1];
sx q[1];
rz(-1.0449303) q[1];
sx q[1];
rz(-2.9108564) q[1];
rz(-pi) q[2];
rz(2.1319904) q[3];
sx q[3];
rz(-1.5586886) q[3];
sx q[3];
rz(-2.5502513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.052496584) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(2.055114) q[2];
rz(2.1975482) q[3];
sx q[3];
rz(-0.090525301) q[3];
sx q[3];
rz(2.0199147) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1679967) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(-2.9214389) q[0];
rz(-1.5036229) q[1];
sx q[1];
rz(-1.817768) q[1];
sx q[1];
rz(-2.5880609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93388825) q[0];
sx q[0];
rz(-1.3739062) q[0];
sx q[0];
rz(-1.6879199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4398063) q[2];
sx q[2];
rz(-0.39798073) q[2];
sx q[2];
rz(-0.27962886) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11343918) q[1];
sx q[1];
rz(-0.69481297) q[1];
sx q[1];
rz(3.0656612) q[1];
rz(-1.7995076) q[3];
sx q[3];
rz(-0.89806496) q[3];
sx q[3];
rz(-1.3474238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0452051) q[2];
sx q[2];
rz(-1.4676981) q[2];
sx q[2];
rz(-2.9333147) q[2];
rz(-1.1194718) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49523062) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(1.0736504) q[0];
rz(-3.011009) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(-2.1317587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9511418) q[0];
sx q[0];
rz(-1.1755687) q[0];
sx q[0];
rz(-1.7717096) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6442579) q[2];
sx q[2];
rz(-1.7188004) q[2];
sx q[2];
rz(2.5876887) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6679935) q[1];
sx q[1];
rz(-0.58071857) q[1];
sx q[1];
rz(0.16814997) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4773682) q[3];
sx q[3];
rz(-1.281732) q[3];
sx q[3];
rz(-2.5632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4449571) q[2];
sx q[2];
rz(-0.37142763) q[2];
sx q[2];
rz(-2.8182287) q[2];
rz(-2.4671386) q[3];
sx q[3];
rz(-2.2489397) q[3];
sx q[3];
rz(3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12482878) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-1.0892518) q[0];
rz(2.543154) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(-2.3425897) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865144) q[0];
sx q[0];
rz(-1.8505972) q[0];
sx q[0];
rz(2.8248019) q[0];
rz(0.094303294) q[2];
sx q[2];
rz(-0.906773) q[2];
sx q[2];
rz(1.2066598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5889783) q[1];
sx q[1];
rz(-1.8726906) q[1];
sx q[1];
rz(-1.4149069) q[1];
rz(2.899762) q[3];
sx q[3];
rz(-1.5655091) q[3];
sx q[3];
rz(-1.2323123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4716855) q[2];
sx q[2];
rz(-2.6354463) q[2];
sx q[2];
rz(1.8890107) q[2];
rz(2.0910828) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(2.2091776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893386) q[0];
sx q[0];
rz(-1.4006571) q[0];
sx q[0];
rz(-0.5156714) q[0];
rz(1.4562666) q[1];
sx q[1];
rz(-1.8003502) q[1];
sx q[1];
rz(2.8175443) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0312248) q[0];
sx q[0];
rz(-2.661781) q[0];
sx q[0];
rz(-0.57439248) q[0];
x q[1];
rz(2.6921047) q[2];
sx q[2];
rz(-2.1725093) q[2];
sx q[2];
rz(-2.4488317) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16897136) q[1];
sx q[1];
rz(-1.4310799) q[1];
sx q[1];
rz(2.0129544) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1284655) q[3];
sx q[3];
rz(-1.6383759) q[3];
sx q[3];
rz(2.9014719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1590165) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(-2.4328361) q[2];
rz(-0.7835663) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(-1.6243352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6056972) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(2.9851483) q[0];
rz(-0.63977301) q[1];
sx q[1];
rz(-2.4867058) q[1];
sx q[1];
rz(-0.99008647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759686) q[0];
sx q[0];
rz(-2.1556615) q[0];
sx q[0];
rz(-2.3192899) q[0];
rz(-0.089870139) q[2];
sx q[2];
rz(-1.9192358) q[2];
sx q[2];
rz(-0.2350829) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53452021) q[1];
sx q[1];
rz(-1.7036408) q[1];
sx q[1];
rz(-0.55664363) q[1];
rz(1.7016497) q[3];
sx q[3];
rz(-1.3706809) q[3];
sx q[3];
rz(-0.36267422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5797552) q[2];
sx q[2];
rz(-1.4116776) q[2];
sx q[2];
rz(-0.61335316) q[2];
rz(-2.8776045) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-0.67155513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494674) q[0];
sx q[0];
rz(-1.6138074) q[0];
sx q[0];
rz(-2.001979) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(1.175066) q[2];
sx q[2];
rz(-1.9446951) q[2];
sx q[2];
rz(-2.4315628) q[2];
rz(-3.0861978) q[3];
sx q[3];
rz(-0.61943409) q[3];
sx q[3];
rz(-2.7243765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
