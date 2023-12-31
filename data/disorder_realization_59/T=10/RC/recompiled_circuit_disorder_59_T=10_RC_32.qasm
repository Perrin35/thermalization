OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(0.27702364) q[1];
sx q[1];
rz(-0.47200534) q[1];
sx q[1];
rz(-3.1402052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35933094) q[0];
sx q[0];
rz(-1.2027272) q[0];
sx q[0];
rz(2.9666535) q[0];
x q[1];
rz(-3.0543047) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(-1.0686312) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5056155) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(-0.14640267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9790768) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(-2.0548267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15443054) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(2.7367676) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(1.0372112) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(2.326139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.095123) q[0];
sx q[0];
rz(-2.2182584) q[0];
sx q[0];
rz(1.4955273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4120861) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(1.6934998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(-1.5335598) q[1];
rz(-pi) q[2];
rz(0.64579441) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(-1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-2.2667623) q[0];
rz(-1.8114999) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(-0.99951807) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(-0.15655984) q[0];
rz(2.2314084) q[2];
sx q[2];
rz(-1.0599469) q[2];
sx q[2];
rz(0.78188932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1516583) q[1];
sx q[1];
rz(-0.5608359) q[1];
sx q[1];
rz(-0.49353091) q[1];
rz(-pi) q[2];
rz(-2.0796892) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(-3.0825465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(-0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(0.27483637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7461473) q[0];
sx q[0];
rz(-1.4822042) q[0];
sx q[0];
rz(-0.079061411) q[0];
rz(-1.5675797) q[2];
sx q[2];
rz(-0.46041691) q[2];
sx q[2];
rz(2.761063) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53510016) q[1];
sx q[1];
rz(-0.98787687) q[1];
sx q[1];
rz(-0.54775723) q[1];
rz(-3.0148388) q[3];
sx q[3];
rz(-0.8781913) q[3];
sx q[3];
rz(2.9664489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(-0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18407962) q[0];
sx q[0];
rz(-0.74680579) q[0];
sx q[0];
rz(2.1225131) q[0];
rz(0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(-1.595572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9285674) q[1];
sx q[1];
rz(-1.2270317) q[1];
sx q[1];
rz(2.7108971) q[1];
rz(-0.77002854) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(0.3240164) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(1.5312622) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76535392) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37939385) q[0];
sx q[0];
rz(-1.6025935) q[0];
sx q[0];
rz(1.6367153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2773877) q[2];
sx q[2];
rz(-0.83653203) q[2];
sx q[2];
rz(-2.4794527) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.016547116) q[1];
sx q[1];
rz(-2.6066337) q[1];
sx q[1];
rz(2.6782481) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42231456) q[3];
sx q[3];
rz(-1.1186244) q[3];
sx q[3];
rz(-2.5715668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1910151) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-2.3699956) q[2];
rz(-2.5937882) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(2.5781393) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(0.34564885) q[0];
rz(-3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-0.46494928) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9879887) q[0];
sx q[0];
rz(-2.0360332) q[0];
sx q[0];
rz(2.9265755) q[0];
x q[1];
rz(-0.49352383) q[2];
sx q[2];
rz(-1.0527305) q[2];
sx q[2];
rz(1.354419) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9629434) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(-0.21957285) q[1];
x q[2];
rz(-2.081359) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(-0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1348304) q[0];
sx q[0];
rz(-1.6400596) q[0];
sx q[0];
rz(0.80471054) q[0];
rz(-pi) q[1];
rz(2.3381091) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(1.4575046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2625418) q[1];
sx q[1];
rz(-1.3971551) q[1];
sx q[1];
rz(0.91211984) q[1];
rz(-pi) q[2];
rz(0.034268495) q[3];
sx q[3];
rz(-1.4141603) q[3];
sx q[3];
rz(-2.6019415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-2.233278) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(0.58319432) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(-0.87402469) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6459991) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(1.5587224) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85272312) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(-0.66327099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45174949) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(-1.4228574) q[1];
x q[2];
rz(-0.60871082) q[3];
sx q[3];
rz(-1.8551644) q[3];
sx q[3];
rz(1.1876719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.727227) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282384) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(0.21270919) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4775425) q[2];
sx q[2];
rz(-1.6972731) q[2];
sx q[2];
rz(-1.7699514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.57943343) q[1];
sx q[1];
rz(-0.65083083) q[1];
sx q[1];
rz(-1.2594373) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3801109) q[3];
sx q[3];
rz(-2.5862525) q[3];
sx q[3];
rz(-1.4243319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5205004) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(2.3836366) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(2.010871) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(1.422613) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
