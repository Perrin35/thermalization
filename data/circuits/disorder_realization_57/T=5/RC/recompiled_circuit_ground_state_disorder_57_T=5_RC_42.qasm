OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(0.33696365) q[0];
rz(-2.872074) q[1];
sx q[1];
rz(-0.94269204) q[1];
sx q[1];
rz(1.8820794) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8963321) q[0];
sx q[0];
rz(-1.0335644) q[0];
sx q[0];
rz(2.2396829) q[0];
rz(-pi) q[1];
rz(-1.7808229) q[2];
sx q[2];
rz(-0.84473306) q[2];
sx q[2];
rz(0.10440102) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14447695) q[1];
sx q[1];
rz(-2.1142748) q[1];
sx q[1];
rz(-1.687084) q[1];
rz(-0.41347031) q[3];
sx q[3];
rz(-1.0116825) q[3];
sx q[3];
rz(-2.632189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9691465) q[2];
sx q[2];
rz(-1.5643876) q[2];
sx q[2];
rz(-1.9077612) q[2];
rz(1.5305758) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(-0.56510258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81561404) q[0];
sx q[0];
rz(-2.1765206) q[0];
sx q[0];
rz(0.24398971) q[0];
rz(-1.1583534) q[1];
sx q[1];
rz(-1.0141677) q[1];
sx q[1];
rz(-2.6932531) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8566003) q[0];
sx q[0];
rz(-2.0584848) q[0];
sx q[0];
rz(1.7708805) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2107606) q[2];
sx q[2];
rz(-2.2907567) q[2];
sx q[2];
rz(-2.623327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92384752) q[1];
sx q[1];
rz(-0.73647803) q[1];
sx q[1];
rz(-2.272241) q[1];
x q[2];
rz(2.118894) q[3];
sx q[3];
rz(-1.992199) q[3];
sx q[3];
rz(0.92589007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1925194) q[2];
sx q[2];
rz(-3.101888) q[2];
sx q[2];
rz(0.81083361) q[2];
rz(-0.0090946322) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(0.31753376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72967616) q[0];
sx q[0];
rz(-1.7971973) q[0];
sx q[0];
rz(-0.94773951) q[0];
rz(-2.2432227) q[1];
sx q[1];
rz(-2.4285474) q[1];
sx q[1];
rz(-0.20634849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2437993) q[0];
sx q[0];
rz(-2.3256105) q[0];
sx q[0];
rz(3.131098) q[0];
rz(1.1694109) q[2];
sx q[2];
rz(-0.81935173) q[2];
sx q[2];
rz(1.2945021) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.628133) q[1];
sx q[1];
rz(-2.0230789) q[1];
sx q[1];
rz(0.22925218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.815229) q[3];
sx q[3];
rz(-1.6152528) q[3];
sx q[3];
rz(1.8281368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72948939) q[2];
sx q[2];
rz(-2.1495843) q[2];
sx q[2];
rz(-0.98423973) q[2];
rz(-0.49736831) q[3];
sx q[3];
rz(-1.2593185) q[3];
sx q[3];
rz(-0.74500144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3618149) q[0];
sx q[0];
rz(-0.091876939) q[0];
sx q[0];
rz(0.20633695) q[0];
rz(1.4641209) q[1];
sx q[1];
rz(-2.3787777) q[1];
sx q[1];
rz(0.62613553) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0968321) q[0];
sx q[0];
rz(-2.2683168) q[0];
sx q[0];
rz(-2.6768854) q[0];
rz(-pi) q[1];
rz(-3.1244833) q[2];
sx q[2];
rz(-0.10280156) q[2];
sx q[2];
rz(-0.67827618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4724222) q[1];
sx q[1];
rz(-2.0991994) q[1];
sx q[1];
rz(-2.1329358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3317165) q[3];
sx q[3];
rz(-2.7446973) q[3];
sx q[3];
rz(-1.797685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76554406) q[2];
sx q[2];
rz(-0.77771336) q[2];
sx q[2];
rz(-0.98975873) q[2];
rz(1.0342213) q[3];
sx q[3];
rz(-1.9725622) q[3];
sx q[3];
rz(-2.6160252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4393716) q[0];
sx q[0];
rz(-1.7218497) q[0];
sx q[0];
rz(-1.9787582) q[0];
rz(-0.16712664) q[1];
sx q[1];
rz(-1.6163328) q[1];
sx q[1];
rz(0.67536813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9769558) q[0];
sx q[0];
rz(-2.2222493) q[0];
sx q[0];
rz(1.8465592) q[0];
rz(-pi) q[1];
rz(-0.33950342) q[2];
sx q[2];
rz(-1.9304747) q[2];
sx q[2];
rz(0.24696697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6847303) q[1];
sx q[1];
rz(-2.5764675) q[1];
sx q[1];
rz(-0.88688382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9642595) q[3];
sx q[3];
rz(-4*pi/5) q[3];
sx q[3];
rz(-2.2817557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9250138) q[2];
sx q[2];
rz(-1.719097) q[2];
sx q[2];
rz(0.50755802) q[2];
rz(3.1223068) q[3];
sx q[3];
rz(-2.3465893) q[3];
sx q[3];
rz(1.9222586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60292202) q[0];
sx q[0];
rz(-2.5141073) q[0];
sx q[0];
rz(2.4399309) q[0];
rz(-2.498846) q[1];
sx q[1];
rz(-1.9822491) q[1];
sx q[1];
rz(1.3210375) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3202527) q[0];
sx q[0];
rz(-0.62766111) q[0];
sx q[0];
rz(1.5492803) q[0];
x q[1];
rz(2.7847544) q[2];
sx q[2];
rz(-2.2388487) q[2];
sx q[2];
rz(2.2989863) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9737742) q[1];
sx q[1];
rz(-1.201448) q[1];
sx q[1];
rz(-1.0008873) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28863971) q[3];
sx q[3];
rz(-0.83191365) q[3];
sx q[3];
rz(-0.4543685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8464437) q[2];
sx q[2];
rz(-0.48131338) q[2];
sx q[2];
rz(0.63977891) q[2];
rz(-0.65711895) q[3];
sx q[3];
rz(-1.6491456) q[3];
sx q[3];
rz(2.8583728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101864) q[0];
sx q[0];
rz(-1.2861847) q[0];
sx q[0];
rz(-0.78045994) q[0];
rz(1.2377493) q[1];
sx q[1];
rz(-1.8433808) q[1];
sx q[1];
rz(2.9775528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1262948) q[0];
sx q[0];
rz(-1.8079213) q[0];
sx q[0];
rz(-2.1046647) q[0];
rz(-pi) q[1];
rz(-0.41167792) q[2];
sx q[2];
rz(-1.9814166) q[2];
sx q[2];
rz(-0.1146929) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8143927) q[1];
sx q[1];
rz(-0.76451028) q[1];
sx q[1];
rz(2.3491377) q[1];
x q[2];
rz(2.3497054) q[3];
sx q[3];
rz(-1.2372176) q[3];
sx q[3];
rz(-1.4871979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26874545) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(2.2609113) q[2];
rz(1.667048) q[3];
sx q[3];
rz(-2.0737952) q[3];
sx q[3];
rz(-1.7269945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47377652) q[0];
sx q[0];
rz(-0.18874636) q[0];
sx q[0];
rz(-1.1291946) q[0];
rz(-0.95343626) q[1];
sx q[1];
rz(-0.9402746) q[1];
sx q[1];
rz(0.58852351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1746219) q[0];
sx q[0];
rz(-2.1123288) q[0];
sx q[0];
rz(-2.7906899) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2382201) q[2];
sx q[2];
rz(-1.7875449) q[2];
sx q[2];
rz(0.34453604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1037081) q[1];
sx q[1];
rz(-2.3988535) q[1];
sx q[1];
rz(0.4813511) q[1];
x q[2];
rz(0.053003691) q[3];
sx q[3];
rz(-1.7946968) q[3];
sx q[3];
rz(-0.25389029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15527655) q[2];
sx q[2];
rz(-2.1337815) q[2];
sx q[2];
rz(1.0922095) q[2];
rz(-2.3327667) q[3];
sx q[3];
rz(-1.0816962) q[3];
sx q[3];
rz(1.697418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5876193) q[0];
sx q[0];
rz(-1.076979) q[0];
sx q[0];
rz(2.4699566) q[0];
rz(2.6486168) q[1];
sx q[1];
rz(-2.2608345) q[1];
sx q[1];
rz(-1.902045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15970358) q[0];
sx q[0];
rz(-0.23724876) q[0];
sx q[0];
rz(1.5818198) q[0];
rz(2.9443594) q[2];
sx q[2];
rz(-1.5610362) q[2];
sx q[2];
rz(0.41744864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18136637) q[1];
sx q[1];
rz(-2.005026) q[1];
sx q[1];
rz(-0.63970345) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9734816) q[3];
sx q[3];
rz(-0.77606499) q[3];
sx q[3];
rz(-0.69918442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.12941831) q[2];
sx q[2];
rz(-1.6678026) q[2];
sx q[2];
rz(-2.5950281) q[2];
rz(0.92285815) q[3];
sx q[3];
rz(-2.512629) q[3];
sx q[3];
rz(-1.1465237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4839812) q[0];
sx q[0];
rz(-0.46009362) q[0];
sx q[0];
rz(2.7125603) q[0];
rz(1.8589164) q[1];
sx q[1];
rz(-1.7762643) q[1];
sx q[1];
rz(-0.081258953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6264296) q[0];
sx q[0];
rz(-1.4019971) q[0];
sx q[0];
rz(-1.6506881) q[0];
rz(-pi) q[1];
rz(-0.0325412) q[2];
sx q[2];
rz(-1.1109118) q[2];
sx q[2];
rz(-2.3655917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4119307) q[1];
sx q[1];
rz(-2.2457128) q[1];
sx q[1];
rz(2.6107671) q[1];
x q[2];
rz(0.83286442) q[3];
sx q[3];
rz(-0.64117764) q[3];
sx q[3];
rz(0.26713757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73605865) q[2];
sx q[2];
rz(-0.93090504) q[2];
sx q[2];
rz(-2.0299358) q[2];
rz(-1.9723655) q[3];
sx q[3];
rz(-2.4284095) q[3];
sx q[3];
rz(-2.9986103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0040759857) q[0];
sx q[0];
rz(-0.62340323) q[0];
sx q[0];
rz(2.2829983) q[0];
rz(-0.89523347) q[1];
sx q[1];
rz(-1.8276855) q[1];
sx q[1];
rz(-3.102416) q[1];
rz(-1.3139541) q[2];
sx q[2];
rz(-1.5594635) q[2];
sx q[2];
rz(2.9008627) q[2];
rz(0.64057401) q[3];
sx q[3];
rz(-1.4645897) q[3];
sx q[3];
rz(0.80932643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
