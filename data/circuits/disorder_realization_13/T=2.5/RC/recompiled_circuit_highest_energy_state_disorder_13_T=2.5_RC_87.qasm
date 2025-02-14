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
rz(0.27999347) q[0];
sx q[0];
rz(-1.0323098) q[0];
sx q[0];
rz(-1.1601467) q[0];
rz(-0.87814826) q[1];
sx q[1];
rz(-2.703009) q[1];
sx q[1];
rz(2.261472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11186799) q[0];
sx q[0];
rz(-0.97767413) q[0];
sx q[0];
rz(0.2082227) q[0];
x q[1];
rz(1.1678668) q[2];
sx q[2];
rz(-1.4968728) q[2];
sx q[2];
rz(-1.5396724) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62978201) q[1];
sx q[1];
rz(-1.7336048) q[1];
sx q[1];
rz(0.96624272) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4515458) q[3];
sx q[3];
rz(-2.0219943) q[3];
sx q[3];
rz(2.638272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1151513) q[2];
sx q[2];
rz(-2.5796964) q[2];
sx q[2];
rz(2.0217516) q[2];
rz(2.1378873) q[3];
sx q[3];
rz(-2.7916838) q[3];
sx q[3];
rz(-2.2899535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.411946) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(-0.44806421) q[0];
rz(-1.0753151) q[1];
sx q[1];
rz(-1.0766462) q[1];
sx q[1];
rz(3.101128) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4498398) q[0];
sx q[0];
rz(-1.7621855) q[0];
sx q[0];
rz(2.2474849) q[0];
rz(1.984111) q[2];
sx q[2];
rz(-1.6530452) q[2];
sx q[2];
rz(-2.8678107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76772122) q[1];
sx q[1];
rz(-1.5295856) q[1];
sx q[1];
rz(-0.26880625) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026626066) q[3];
sx q[3];
rz(-2.8288719) q[3];
sx q[3];
rz(-1.3440162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2761817) q[2];
sx q[2];
rz(-0.9006331) q[2];
sx q[2];
rz(0.5033699) q[2];
rz(1.447575) q[3];
sx q[3];
rz(-1.9083128) q[3];
sx q[3];
rz(2.2342822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4564948) q[0];
sx q[0];
rz(-2.1879897) q[0];
sx q[0];
rz(-2.8431235) q[0];
rz(-1.5553156) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(-1.5063937) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54558447) q[0];
sx q[0];
rz(-2.8337361) q[0];
sx q[0];
rz(2.8915384) q[0];
rz(-pi) q[1];
rz(0.2770642) q[2];
sx q[2];
rz(-2.2383406) q[2];
sx q[2];
rz(0.88457047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7506444) q[1];
sx q[1];
rz(-1.4381289) q[1];
sx q[1];
rz(0.68718095) q[1];
rz(-pi) q[2];
rz(-1.8465145) q[3];
sx q[3];
rz(-2.9447768) q[3];
sx q[3];
rz(-2.4829602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35206587) q[2];
sx q[2];
rz(-2.3730998) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(-1.7302892) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(2.3119149) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4378433) q[0];
sx q[0];
rz(-2.4898536) q[0];
sx q[0];
rz(-1.334345) q[0];
rz(-1.5127381) q[1];
sx q[1];
rz(-1.6213497) q[1];
sx q[1];
rz(-0.13660647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2481987) q[0];
sx q[0];
rz(-1.8780439) q[0];
sx q[0];
rz(0.91715468) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2567472) q[2];
sx q[2];
rz(-2.8171091) q[2];
sx q[2];
rz(2.9195234) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1699511) q[1];
sx q[1];
rz(-1.2593237) q[1];
sx q[1];
rz(0.71572742) q[1];
rz(-pi) q[2];
rz(2.8084561) q[3];
sx q[3];
rz(-1.9018947) q[3];
sx q[3];
rz(-3.0931728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85285774) q[2];
sx q[2];
rz(-1.6625657) q[2];
sx q[2];
rz(-0.23957254) q[2];
rz(-2.1824956) q[3];
sx q[3];
rz(-1.164091) q[3];
sx q[3];
rz(-1.0378999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94944823) q[0];
sx q[0];
rz(-1.0553772) q[0];
sx q[0];
rz(-1.5437833) q[0];
rz(-1.1100618) q[1];
sx q[1];
rz(-1.9381356) q[1];
sx q[1];
rz(-1.69453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2758165) q[0];
sx q[0];
rz(-1.5636684) q[0];
sx q[0];
rz(0.29840666) q[0];
x q[1];
rz(1.6844774) q[2];
sx q[2];
rz(-0.34514134) q[2];
sx q[2];
rz(-0.98649401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.49741433) q[1];
sx q[1];
rz(-1.7458054) q[1];
sx q[1];
rz(2.464556) q[1];
rz(1.9915699) q[3];
sx q[3];
rz(-2.07486) q[3];
sx q[3];
rz(-2.1114388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99628228) q[2];
sx q[2];
rz(-2.0029533) q[2];
sx q[2];
rz(-0.81926695) q[2];
rz(0.94775689) q[3];
sx q[3];
rz(-2.3534333) q[3];
sx q[3];
rz(-1.88571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6325697) q[0];
sx q[0];
rz(-3.0036354) q[0];
sx q[0];
rz(-2.5914958) q[0];
rz(0.24712786) q[1];
sx q[1];
rz(-1.1526266) q[1];
sx q[1];
rz(-2.1955042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23299117) q[0];
sx q[0];
rz(-2.6938022) q[0];
sx q[0];
rz(-0.93771387) q[0];
rz(-pi) q[1];
rz(-1.2900903) q[2];
sx q[2];
rz(-2.9200313) q[2];
sx q[2];
rz(-2.6803379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4934686) q[1];
sx q[1];
rz(-1.5350684) q[1];
sx q[1];
rz(-2.6717527) q[1];
rz(-pi) q[2];
rz(0.23885085) q[3];
sx q[3];
rz(-1.672604) q[3];
sx q[3];
rz(0.073697173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3863824) q[2];
sx q[2];
rz(-2.5911665) q[2];
sx q[2];
rz(0.68106252) q[2];
rz(0.81124535) q[3];
sx q[3];
rz(-2.097605) q[3];
sx q[3];
rz(0.55027858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193467) q[0];
sx q[0];
rz(-0.62048727) q[0];
sx q[0];
rz(0.74561179) q[0];
rz(1.1634119) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(-2.9643639) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39800571) q[0];
sx q[0];
rz(-1.2400125) q[0];
sx q[0];
rz(2.8547006) q[0];
rz(-pi) q[1];
rz(1.6615599) q[2];
sx q[2];
rz(-1.5077733) q[2];
sx q[2];
rz(-1.3985967) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2043269) q[1];
sx q[1];
rz(-2.4961053) q[1];
sx q[1];
rz(0.93475229) q[1];
rz(-0.69474649) q[3];
sx q[3];
rz(-1.5383895) q[3];
sx q[3];
rz(2.29914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54997286) q[2];
sx q[2];
rz(-2.1592996) q[2];
sx q[2];
rz(1.4493235) q[2];
rz(2.5189404) q[3];
sx q[3];
rz(-2.6144274) q[3];
sx q[3];
rz(-1.8221633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016668884) q[0];
sx q[0];
rz(-1.7542087) q[0];
sx q[0];
rz(1.5119875) q[0];
rz(-2.2177057) q[1];
sx q[1];
rz(-1.7216564) q[1];
sx q[1];
rz(-2.5249544) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1204471) q[0];
sx q[0];
rz(-2.5207656) q[0];
sx q[0];
rz(1.3070941) q[0];
rz(-pi) q[1];
rz(3.055279) q[2];
sx q[2];
rz(-2.6722135) q[2];
sx q[2];
rz(2.1522107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60395998) q[1];
sx q[1];
rz(-0.62509552) q[1];
sx q[1];
rz(-0.14664607) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8455335) q[3];
sx q[3];
rz(-1.0062075) q[3];
sx q[3];
rz(-0.72455154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7743249) q[2];
sx q[2];
rz(-0.35949817) q[2];
sx q[2];
rz(-2.9765339) q[2];
rz(-0.74602357) q[3];
sx q[3];
rz(-1.6227928) q[3];
sx q[3];
rz(-0.042796854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54166722) q[0];
sx q[0];
rz(-0.2921108) q[0];
sx q[0];
rz(-0.41411972) q[0];
rz(-2.7060624) q[1];
sx q[1];
rz(-0.51594096) q[1];
sx q[1];
rz(1.863265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.965344) q[0];
sx q[0];
rz(-2.8581736) q[0];
sx q[0];
rz(2.118628) q[0];
rz(0.31164557) q[2];
sx q[2];
rz(-1.4053182) q[2];
sx q[2];
rz(-0.32725016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32132253) q[1];
sx q[1];
rz(-1.6386809) q[1];
sx q[1];
rz(1.5589139) q[1];
rz(-pi) q[2];
rz(1.0429562) q[3];
sx q[3];
rz(-1.5164638) q[3];
sx q[3];
rz(-2.899037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4182959) q[2];
sx q[2];
rz(-1.3822684) q[2];
sx q[2];
rz(2.3988775) q[2];
rz(-2.3501979) q[3];
sx q[3];
rz(-2.4580038) q[3];
sx q[3];
rz(-1.4092103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.047569711) q[0];
sx q[0];
rz(-0.3293193) q[0];
sx q[0];
rz(-0.91162115) q[0];
rz(2.1564663) q[1];
sx q[1];
rz(-0.42675012) q[1];
sx q[1];
rz(-1.3381348) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.074756) q[0];
sx q[0];
rz(-1.8843411) q[0];
sx q[0];
rz(-1.0694802) q[0];
rz(-2.2608032) q[2];
sx q[2];
rz(-2.6652209) q[2];
sx q[2];
rz(2.7169189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2240026) q[1];
sx q[1];
rz(-0.75079599) q[1];
sx q[1];
rz(0.96377967) q[1];
rz(-pi) q[2];
rz(2.2959241) q[3];
sx q[3];
rz(-0.4825497) q[3];
sx q[3];
rz(3.1180814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35475874) q[2];
sx q[2];
rz(-1.6666731) q[2];
sx q[2];
rz(1.1653398) q[2];
rz(1.3960086) q[3];
sx q[3];
rz(-1.1895807) q[3];
sx q[3];
rz(-1.5918119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5042481) q[0];
sx q[0];
rz(-1.5550384) q[0];
sx q[0];
rz(0.73873781) q[0];
rz(0.63738102) q[1];
sx q[1];
rz(-1.6045872) q[1];
sx q[1];
rz(2.3789023) q[1];
rz(1.8861072) q[2];
sx q[2];
rz(-2.6985395) q[2];
sx q[2];
rz(0.48816608) q[2];
rz(-0.55841726) q[3];
sx q[3];
rz(-1.1815841) q[3];
sx q[3];
rz(-3.0403137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
