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
rz(2.4861205) q[0];
sx q[0];
rz(-0.95215005) q[0];
sx q[0];
rz(-0.093753554) q[0];
rz(1.3682415) q[1];
sx q[1];
rz(4.7028766) q[1];
sx q[1];
rz(8.75755) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038644636) q[0];
sx q[0];
rz(-1.9460475) q[0];
sx q[0];
rz(-2.1041591) q[0];
x q[1];
rz(2.8831161) q[2];
sx q[2];
rz(-1.7670858) q[2];
sx q[2];
rz(-0.95610134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2161795) q[1];
sx q[1];
rz(-0.88757463) q[1];
sx q[1];
rz(-0.36551468) q[1];
rz(-0.85140075) q[3];
sx q[3];
rz(-2.3150597) q[3];
sx q[3];
rz(-1.9498619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5253456) q[2];
sx q[2];
rz(-1.5364001) q[2];
sx q[2];
rz(-3.1180535) q[2];
rz(-0.25447887) q[3];
sx q[3];
rz(-1.1602217) q[3];
sx q[3];
rz(-2.3833073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94754058) q[0];
sx q[0];
rz(-1.7451311) q[0];
sx q[0];
rz(2.5260455) q[0];
rz(-0.60000348) q[1];
sx q[1];
rz(-1.2702076) q[1];
sx q[1];
rz(2.792865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1407326) q[0];
sx q[0];
rz(-1.494074) q[0];
sx q[0];
rz(2.4984635) q[0];
x q[1];
rz(1.9127513) q[2];
sx q[2];
rz(-1.8445024) q[2];
sx q[2];
rz(-0.051816377) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2243721) q[1];
sx q[1];
rz(-1.1172235) q[1];
sx q[1];
rz(2.808203) q[1];
rz(-1.0660961) q[3];
sx q[3];
rz(-1.821992) q[3];
sx q[3];
rz(-0.78758729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6500924) q[2];
sx q[2];
rz(-1.8365752) q[2];
sx q[2];
rz(0.65518641) q[2];
rz(0.16264597) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(-2.9338525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1732037) q[0];
sx q[0];
rz(-1.116331) q[0];
sx q[0];
rz(0.094245687) q[0];
rz(-1.3000129) q[1];
sx q[1];
rz(-2.2424707) q[1];
sx q[1];
rz(-1.947044) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904063) q[0];
sx q[0];
rz(-1.8739432) q[0];
sx q[0];
rz(2.6180352) q[0];
rz(-pi) q[1];
rz(-0.051812096) q[2];
sx q[2];
rz(-1.9345043) q[2];
sx q[2];
rz(1.6339982) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8988674) q[1];
sx q[1];
rz(-2.8812802) q[1];
sx q[1];
rz(-1.3993466) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1949092) q[3];
sx q[3];
rz(-0.84400153) q[3];
sx q[3];
rz(0.90968859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3208348) q[2];
sx q[2];
rz(-1.4471549) q[2];
sx q[2];
rz(2.7336332) q[2];
rz(1.3784493) q[3];
sx q[3];
rz(-2.2429376) q[3];
sx q[3];
rz(3.0234226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5648062) q[0];
sx q[0];
rz(-1.7498359) q[0];
sx q[0];
rz(0.80192649) q[0];
rz(2.7469514) q[1];
sx q[1];
rz(-2.092974) q[1];
sx q[1];
rz(-3.1065497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3010522) q[0];
sx q[0];
rz(-0.6834712) q[0];
sx q[0];
rz(-2.1816064) q[0];
rz(-1.3801244) q[2];
sx q[2];
rz(-1.241127) q[2];
sx q[2];
rz(-0.53527385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1890244) q[1];
sx q[1];
rz(-0.060563001) q[1];
sx q[1];
rz(2.3409055) q[1];
rz(2.7862042) q[3];
sx q[3];
rz(-1.3878763) q[3];
sx q[3];
rz(-1.2661883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0617712) q[2];
sx q[2];
rz(-2.0786736) q[2];
sx q[2];
rz(-1.4351832) q[2];
rz(1.4052514) q[3];
sx q[3];
rz(-2.0515714) q[3];
sx q[3];
rz(-1.1100356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6768796) q[0];
sx q[0];
rz(-1.9929303) q[0];
sx q[0];
rz(-1.9997464) q[0];
rz(-0.51072085) q[1];
sx q[1];
rz(-1.2341713) q[1];
sx q[1];
rz(1.7150735) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0066442) q[0];
sx q[0];
rz(-2.0922959) q[0];
sx q[0];
rz(-2.5000076) q[0];
rz(-1.2586589) q[2];
sx q[2];
rz(-1.7870047) q[2];
sx q[2];
rz(-2.5463111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74523679) q[1];
sx q[1];
rz(-2.5380683) q[1];
sx q[1];
rz(-1.7698494) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9782039) q[3];
sx q[3];
rz(-1.4168882) q[3];
sx q[3];
rz(-0.28095442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4003754) q[2];
sx q[2];
rz(-1.8343238) q[2];
sx q[2];
rz(0.2571787) q[2];
rz(2.2687965) q[3];
sx q[3];
rz(-1.8200487) q[3];
sx q[3];
rz(2.0040373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0208825) q[0];
sx q[0];
rz(-2.6873984) q[0];
sx q[0];
rz(2.0632451) q[0];
rz(2.2348166) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(-0.041898601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3495765) q[0];
sx q[0];
rz(-2.0431762) q[0];
sx q[0];
rz(0.040332009) q[0];
rz(-2.2380377) q[2];
sx q[2];
rz(-1.4741401) q[2];
sx q[2];
rz(2.5943611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4387233) q[1];
sx q[1];
rz(-1.8449191) q[1];
sx q[1];
rz(-1.1945748) q[1];
rz(-pi) q[2];
rz(-0.53421212) q[3];
sx q[3];
rz(-0.90972661) q[3];
sx q[3];
rz(2.7345757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8580253) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(-0.89135998) q[2];
rz(-2.4274965) q[3];
sx q[3];
rz(-2.4024506) q[3];
sx q[3];
rz(-1.063063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080268) q[0];
sx q[0];
rz(-1.2440246) q[0];
sx q[0];
rz(2.4666393) q[0];
rz(1.630111) q[1];
sx q[1];
rz(-2.5210896) q[1];
sx q[1];
rz(-2.564548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8682315) q[0];
sx q[0];
rz(-1.8087862) q[0];
sx q[0];
rz(2.0711876) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7084025) q[2];
sx q[2];
rz(-1.6640888) q[2];
sx q[2];
rz(-0.46822883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.31957175) q[1];
sx q[1];
rz(-2.4720925) q[1];
sx q[1];
rz(1.7029525) q[1];
rz(-0.8612154) q[3];
sx q[3];
rz(-0.83038721) q[3];
sx q[3];
rz(-1.9423241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6881037) q[2];
sx q[2];
rz(-1.1223531) q[2];
sx q[2];
rz(-2.2171059) q[2];
rz(-1.9214572) q[3];
sx q[3];
rz(-2.852738) q[3];
sx q[3];
rz(-3.0815304) q[3];
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
rz(2.7500551) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(-1.7975988) q[0];
rz(-1.9186107) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(-1.5498243) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90541461) q[0];
sx q[0];
rz(-1.2981725) q[0];
sx q[0];
rz(2.9546066) q[0];
rz(-pi) q[1];
rz(0.28002589) q[2];
sx q[2];
rz(-0.38869263) q[2];
sx q[2];
rz(-0.9946781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0613411) q[1];
sx q[1];
rz(-2.3480573) q[1];
sx q[1];
rz(0.52498753) q[1];
rz(-pi) q[2];
rz(2.377691) q[3];
sx q[3];
rz(-1.502486) q[3];
sx q[3];
rz(2.1602283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0595155) q[2];
sx q[2];
rz(-2.1166708) q[2];
sx q[2];
rz(-2.9151741) q[2];
rz(3.1021127) q[3];
sx q[3];
rz(-1.4791146) q[3];
sx q[3];
rz(-1.1994908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3969642) q[0];
sx q[0];
rz(-1.2122943) q[0];
sx q[0];
rz(-0.82897559) q[0];
rz(-0.12116155) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(-1.4519579) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4273116) q[0];
sx q[0];
rz(-1.7085008) q[0];
sx q[0];
rz(2.4192823) q[0];
rz(0.42134085) q[2];
sx q[2];
rz(-2.2724814) q[2];
sx q[2];
rz(2.7340467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.02579298) q[1];
sx q[1];
rz(-0.83459548) q[1];
sx q[1];
rz(-2.2195312) q[1];
x q[2];
rz(-2.1249346) q[3];
sx q[3];
rz(-1.2264369) q[3];
sx q[3];
rz(0.90079067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5539603) q[2];
sx q[2];
rz(-0.82531896) q[2];
sx q[2];
rz(-2.2770605) q[2];
rz(-0.65166059) q[3];
sx q[3];
rz(-2.103002) q[3];
sx q[3];
rz(3.1237176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5038576) q[0];
sx q[0];
rz(-0.88541579) q[0];
sx q[0];
rz(0.46022415) q[0];
rz(-1.7962615) q[1];
sx q[1];
rz(-2.5049152) q[1];
sx q[1];
rz(1.3705137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11281989) q[0];
sx q[0];
rz(-1.9589073) q[0];
sx q[0];
rz(-0.32916268) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37563328) q[2];
sx q[2];
rz(-0.67258316) q[2];
sx q[2];
rz(-0.022136839) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7975841) q[1];
sx q[1];
rz(-1.3329778) q[1];
sx q[1];
rz(1.306755) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0881365) q[3];
sx q[3];
rz(-1.7776907) q[3];
sx q[3];
rz(-0.89689909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3687849) q[2];
sx q[2];
rz(-1.3489172) q[2];
sx q[2];
rz(0.14300145) q[2];
rz(-0.16608206) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(-0.79296976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.3731257) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(-1.985818) q[1];
sx q[1];
rz(-1.118569) q[1];
sx q[1];
rz(-2.7973693) q[1];
rz(0.58123238) q[2];
sx q[2];
rz(-2.317133) q[2];
sx q[2];
rz(1.9036507) q[2];
rz(1.0977911) q[3];
sx q[3];
rz(-0.65752959) q[3];
sx q[3];
rz(0.65617954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
