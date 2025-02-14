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
rz(-1.8921312) q[0];
rz(3.0013822) q[1];
sx q[1];
rz(-1.4445211) q[1];
sx q[1];
rz(3.0797449) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1309146) q[0];
sx q[0];
rz(-0.82404165) q[0];
sx q[0];
rz(-2.2758621) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.090422618) q[2];
sx q[2];
rz(-1.5704463) q[2];
sx q[2];
rz(-1.4594913) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1143519) q[1];
sx q[1];
rz(-1.4732962) q[1];
sx q[1];
rz(-0.84836324) q[1];
rz(2.2592945) q[3];
sx q[3];
rz(-1.8940482) q[3];
sx q[3];
rz(2.0417449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7141815) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(0.25644914) q[2];
rz(-2.9668258) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8137708) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(-2.9028614) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(3.0366268) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3070827) q[0];
sx q[0];
rz(-1.420546) q[0];
sx q[0];
rz(-0.86850496) q[0];
x q[1];
rz(-0.41609515) q[2];
sx q[2];
rz(-1.7836092) q[2];
sx q[2];
rz(-1.0989604) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2626896) q[1];
sx q[1];
rz(-0.57003747) q[1];
sx q[1];
rz(-1.0703342) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41173287) q[3];
sx q[3];
rz(-0.74854088) q[3];
sx q[3];
rz(2.4553776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77439848) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(1.6801838) q[2];
rz(1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(-0.27209601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-0.21037924) q[0];
sx q[0];
rz(-2.126597) q[0];
rz(-0.89730942) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(0.6764594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7583026) q[0];
sx q[0];
rz(-1.5448017) q[0];
sx q[0];
rz(-2.1878408) q[0];
rz(0.93134201) q[2];
sx q[2];
rz(-0.73340511) q[2];
sx q[2];
rz(-0.99898192) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30005042) q[1];
sx q[1];
rz(-1.9498697) q[1];
sx q[1];
rz(-2.9247523) q[1];
rz(-pi) q[2];
rz(0.81636103) q[3];
sx q[3];
rz(-1.584017) q[3];
sx q[3];
rz(-0.17133443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2066388) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(0.76033956) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-0.81086603) q[3];
sx q[3];
rz(-0.84347239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.7059785) q[0];
sx q[0];
rz(-2.5137081) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(-0.71890038) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(-0.23153201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0569772) q[0];
sx q[0];
rz(-1.9997134) q[0];
sx q[0];
rz(-0.57022988) q[0];
rz(0.79265742) q[2];
sx q[2];
rz(-1.5680299) q[2];
sx q[2];
rz(-2.1994928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0147103) q[1];
sx q[1];
rz(-2.0538804) q[1];
sx q[1];
rz(1.6761858) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9636964) q[3];
sx q[3];
rz(-2.1390954) q[3];
sx q[3];
rz(-1.6015805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71451688) q[2];
sx q[2];
rz(-0.38893739) q[2];
sx q[2];
rz(-0.60849774) q[2];
rz(0.19615873) q[3];
sx q[3];
rz(-1.720263) q[3];
sx q[3];
rz(0.99854809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683559) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-2.3625145) q[0];
rz(2.211606) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.1735865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1176276) q[0];
sx q[0];
rz(-1.574734) q[0];
sx q[0];
rz(1.9088452) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1861325) q[2];
sx q[2];
rz(-1.460607) q[2];
sx q[2];
rz(2.2249391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1184349) q[1];
sx q[1];
rz(-1.6350606) q[1];
sx q[1];
rz(-1.2530909) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7608123) q[3];
sx q[3];
rz(-1.992986) q[3];
sx q[3];
rz(-1.208153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8194627) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(-0.67506153) q[2];
rz(0.86554646) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(1.2077829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753767) q[0];
sx q[0];
rz(-3.1238811) q[0];
sx q[0];
rz(2.2702763) q[0];
rz(-0.21698347) q[1];
sx q[1];
rz(-1.5996409) q[1];
sx q[1];
rz(1.3380231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8080475) q[0];
sx q[0];
rz(-1.6535954) q[0];
sx q[0];
rz(-2.4890578) q[0];
x q[1];
rz(0.14071847) q[2];
sx q[2];
rz(-2.291516) q[2];
sx q[2];
rz(0.36397935) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1068145) q[1];
sx q[1];
rz(-0.15131525) q[1];
sx q[1];
rz(-1.8968015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35049482) q[3];
sx q[3];
rz(-0.51096254) q[3];
sx q[3];
rz(3.0929589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.013082144) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(2.3415372) q[2];
rz(-2.1498146) q[3];
sx q[3];
rz(-2.1938775) q[3];
sx q[3];
rz(-2.1984524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8868788) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(-0.15175858) q[0];
rz(-1.4166547) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-0.83622611) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9660824) q[0];
sx q[0];
rz(-1.9608736) q[0];
sx q[0];
rz(-0.83781029) q[0];
rz(-pi) q[1];
rz(-1.56244) q[2];
sx q[2];
rz(-0.99412336) q[2];
sx q[2];
rz(-1.6461262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9735896) q[1];
sx q[1];
rz(-2.8781946) q[1];
sx q[1];
rz(1.9067287) q[1];
rz(-pi) q[2];
rz(1.4400021) q[3];
sx q[3];
rz(-2.0328641) q[3];
sx q[3];
rz(2.3123284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4862711) q[2];
sx q[2];
rz(-3.0781015) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(2.4231353) q[3];
sx q[3];
rz(-1.5114762) q[3];
sx q[3];
rz(0.92459905) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543168) q[0];
sx q[0];
rz(-0.59816718) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(-0.68880853) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(0.93592962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0416695) q[0];
sx q[0];
rz(-1.3893621) q[0];
sx q[0];
rz(2.6122841) q[0];
x q[1];
rz(1.9804968) q[2];
sx q[2];
rz(-1.8852703) q[2];
sx q[2];
rz(-1.0057698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3080146) q[1];
sx q[1];
rz(-1.2460099) q[1];
sx q[1];
rz(-3.0851425) q[1];
rz(-pi) q[2];
rz(-0.0015663107) q[3];
sx q[3];
rz(-0.74347444) q[3];
sx q[3];
rz(1.8494128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9800637) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(-2.8847983) q[2];
rz(-2.1234546) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7629906) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(2.2879404) q[0];
rz(-1.8866106) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(-2.8010211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2947663) q[0];
sx q[0];
rz(-1.7309523) q[0];
sx q[0];
rz(1.7677884) q[0];
rz(1.3236295) q[2];
sx q[2];
rz(-0.91193141) q[2];
sx q[2];
rz(2.5160922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7813959) q[1];
sx q[1];
rz(-2.2193647) q[1];
sx q[1];
rz(-0.99332033) q[1];
rz(-1.1084308) q[3];
sx q[3];
rz(-1.8183072) q[3];
sx q[3];
rz(0.68456291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7178932) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(0.37810668) q[2];
rz(-0.2374436) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(0.28089359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963999) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(-2.4943446) q[0];
rz(-1.9735533) q[1];
sx q[1];
rz(-0.85473514) q[1];
sx q[1];
rz(-0.38356575) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952571) q[0];
sx q[0];
rz(-1.5492166) q[0];
sx q[0];
rz(3.1109592) q[0];
x q[1];
rz(2.6099989) q[2];
sx q[2];
rz(-1.8583014) q[2];
sx q[2];
rz(0.93898857) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3473222) q[1];
sx q[1];
rz(-1.8350661) q[1];
sx q[1];
rz(-1.3572973) q[1];
rz(-pi) q[2];
rz(1.7586772) q[3];
sx q[3];
rz(-2.6184824) q[3];
sx q[3];
rz(-2.1818102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26455227) q[2];
sx q[2];
rz(-0.88901797) q[2];
sx q[2];
rz(1.5846579) q[2];
rz(-1.5133739) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-0.42678601) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.1813286) q[1];
sx q[1];
rz(-1.7971296) q[1];
sx q[1];
rz(2.8167579) q[1];
rz(0.2427553) q[2];
sx q[2];
rz(-1.1052255) q[2];
sx q[2];
rz(0.33585264) q[2];
rz(-3.0814617) q[3];
sx q[3];
rz(-0.48532607) q[3];
sx q[3];
rz(2.8564712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
