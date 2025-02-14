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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(1.2494614) q[0];
rz(-3.2818031) q[1];
sx q[1];
rz(1.6970716) q[1];
sx q[1];
rz(12.628218) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1309146) q[0];
sx q[0];
rz(-0.82404165) q[0];
sx q[0];
rz(2.2758621) q[0];
rz(-3.1377162) q[2];
sx q[2];
rz(-0.090423294) q[2];
sx q[2];
rz(-0.10744444) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.34644352) q[1];
sx q[1];
rz(-0.7278053) q[1];
sx q[1];
rz(1.4239356) q[1];
x q[2];
rz(-2.7322686) q[3];
sx q[3];
rz(-0.92420368) q[3];
sx q[3];
rz(-2.4149946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7141815) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(-2.8851435) q[2];
rz(0.17476684) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(2.9028614) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(0.10496584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83451) q[0];
sx q[0];
rz(-1.420546) q[0];
sx q[0];
rz(0.86850496) q[0];
rz(0.41609515) q[2];
sx q[2];
rz(-1.3579835) q[2];
sx q[2];
rz(-1.0989604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30279494) q[1];
sx q[1];
rz(-1.0775591) q[1];
sx q[1];
rz(-0.29839813) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2148923) q[3];
sx q[3];
rz(-0.89733636) q[3];
sx q[3];
rz(1.9178903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3671942) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(1.6801838) q[2];
rz(1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(2.126597) q[0];
rz(-0.89730942) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(-0.6764594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7583026) q[0];
sx q[0];
rz(-1.5448017) q[0];
sx q[0];
rz(0.95375188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6482267) q[2];
sx q[2];
rz(-1.0037494) q[2];
sx q[2];
rz(2.9285448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8415422) q[1];
sx q[1];
rz(-1.9498697) q[1];
sx q[1];
rz(-2.9247523) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1234497) q[3];
sx q[3];
rz(-2.3251495) q[3];
sx q[3];
rz(1.7545561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93495381) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(0.76033956) q[2];
rz(1.2065411) q[3];
sx q[3];
rz(-0.81086603) q[3];
sx q[3];
rz(0.84347239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059785) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(0.71890038) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(2.9100606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0615094) q[0];
sx q[0];
rz(-2.4426547) q[0];
sx q[0];
rz(2.438758) q[0];
rz(-pi) q[1];
rz(1.5747373) q[2];
sx q[2];
rz(-0.77814279) q[2];
sx q[2];
rz(-2.5100894) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0147103) q[1];
sx q[1];
rz(-1.0877123) q[1];
sx q[1];
rz(1.6761858) q[1];
rz(-pi) q[2];
rz(-0.54010424) q[3];
sx q[3];
rz(-2.4632235) q[3];
sx q[3];
rz(2.2578007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4270758) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(-0.60849774) q[2];
rz(2.9454339) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(-2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2732368) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-0.77907816) q[0];
rz(-2.211606) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.9680061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0239651) q[0];
sx q[0];
rz(-1.5668586) q[0];
sx q[0];
rz(1.2327475) q[0];
x q[1];
rz(-1.3814244) q[2];
sx q[2];
rz(-0.62386419) q[2];
sx q[2];
rz(-2.3331235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8820937) q[1];
sx q[1];
rz(-2.8176687) q[1];
sx q[1];
rz(1.367635) q[1];
rz(1.7608123) q[3];
sx q[3];
rz(-1.992986) q[3];
sx q[3];
rz(-1.9334396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.32213) q[2];
sx q[2];
rz(-0.16699114) q[2];
sx q[2];
rz(-2.4665311) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38782564) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(2.2702763) q[0];
rz(2.9246092) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(1.8035696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8080475) q[0];
sx q[0];
rz(-1.6535954) q[0];
sx q[0];
rz(-0.6525349) q[0];
rz(-3.0008742) q[2];
sx q[2];
rz(-2.291516) q[2];
sx q[2];
rz(0.36397935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1068145) q[1];
sx q[1];
rz(-0.15131525) q[1];
sx q[1];
rz(1.2447912) q[1];
x q[2];
rz(-2.7910978) q[3];
sx q[3];
rz(-2.6306301) q[3];
sx q[3];
rz(3.0929589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(-0.80005542) q[2];
rz(0.99177805) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(-0.94314027) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8868788) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(-2.9898341) q[0];
rz(1.7249379) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-0.83622611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0043512303) q[0];
sx q[0];
rz(-2.3286331) q[0];
sx q[0];
rz(-2.1217968) q[0];
rz(1.56244) q[2];
sx q[2];
rz(-0.99412336) q[2];
sx q[2];
rz(1.6461262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6265833) q[1];
sx q[1];
rz(-1.8191511) q[1];
sx q[1];
rz(0.088661389) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6760929) q[3];
sx q[3];
rz(-1.6878078) q[3];
sx q[3];
rz(-0.68295628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4862711) q[2];
sx q[2];
rz(-3.0781015) q[2];
sx q[2];
rz(-3.0625694) q[2];
rz(-2.4231353) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(0.92459905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.786161) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(-2.4527841) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(-0.93592962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0416695) q[0];
sx q[0];
rz(-1.7522305) q[0];
sx q[0];
rz(2.6122841) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8008119) q[2];
sx q[2];
rz(-1.182297) q[2];
sx q[2];
rz(-0.69855554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3080146) q[1];
sx q[1];
rz(-1.2460099) q[1];
sx q[1];
rz(-3.0851425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5693561) q[3];
sx q[3];
rz(-2.3142696) q[3];
sx q[3];
rz(-1.8515406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9800637) q[2];
sx q[2];
rz(-2.4825725) q[2];
sx q[2];
rz(-2.8847983) q[2];
rz(2.1234546) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(-2.170678) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37860206) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(-0.85365224) q[0];
rz(1.8866106) q[1];
sx q[1];
rz(-1.9957142) q[1];
sx q[1];
rz(-2.8010211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3981235) q[0];
sx q[0];
rz(-2.8883683) q[0];
sx q[0];
rz(2.2608093) q[0];
rz(-pi) q[1];
rz(-1.8179632) q[2];
sx q[2];
rz(-2.2296612) q[2];
sx q[2];
rz(0.62550046) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53696991) q[1];
sx q[1];
rz(-2.3019362) q[1];
sx q[1];
rz(0.62418681) q[1];
rz(-2.0331618) q[3];
sx q[3];
rz(-1.8183072) q[3];
sx q[3];
rz(2.4570297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7178932) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(-2.763486) q[2];
rz(0.2374436) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(-2.8606991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963999) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(0.64724809) q[0];
rz(1.9735533) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(2.7580269) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18902168) q[0];
sx q[0];
rz(-1.5492166) q[0];
sx q[0];
rz(-0.030633472) q[0];
rz(1.9012544) q[2];
sx q[2];
rz(-1.0631732) q[2];
sx q[2];
rz(-2.3445545) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4870461) q[1];
sx q[1];
rz(-0.33815835) q[1];
sx q[1];
rz(-2.4772993) q[1];
rz(1.0553352) q[3];
sx q[3];
rz(-1.4773507) q[3];
sx q[3];
rz(2.6938113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8770404) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(-1.5569347) q[2];
rz(-1.5133739) q[3];
sx q[3];
rz(-1.3686562) q[3];
sx q[3];
rz(-2.7148066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7814816) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.9602641) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(-2.8988373) q[2];
sx q[2];
rz(-1.1052255) q[2];
sx q[2];
rz(0.33585264) q[2];
rz(-1.539113) q[3];
sx q[3];
rz(-2.0551695) q[3];
sx q[3];
rz(2.7885128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
