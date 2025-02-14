OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1407815) q[0];
sx q[0];
rz(-2.3164764) q[0];
sx q[0];
rz(2.7197279) q[0];
rz(-0.74692625) q[1];
sx q[1];
rz(-2.5513625) q[1];
sx q[1];
rz(2.3205369) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29971545) q[0];
sx q[0];
rz(-1.4869191) q[0];
sx q[0];
rz(1.5257349) q[0];
rz(-pi) q[1];
rz(1.3652322) q[2];
sx q[2];
rz(-1.3495619) q[2];
sx q[2];
rz(1.6541755) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6118879) q[1];
sx q[1];
rz(-2.6571353) q[1];
sx q[1];
rz(0.92682181) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2116634) q[3];
sx q[3];
rz(-2.7756423) q[3];
sx q[3];
rz(-0.32420483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39516285) q[2];
sx q[2];
rz(-2.2300827) q[2];
sx q[2];
rz(2.5968623) q[2];
rz(-1.644545) q[3];
sx q[3];
rz(-2.0293197) q[3];
sx q[3];
rz(-1.3397608) q[3];
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
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(0.5734545) q[0];
rz(0.042512976) q[1];
sx q[1];
rz(-1.7492234) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95623446) q[0];
sx q[0];
rz(-1.6969674) q[0];
sx q[0];
rz(1.5500359) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3147215) q[2];
sx q[2];
rz(-1.8559858) q[2];
sx q[2];
rz(-1.3392753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1309218) q[1];
sx q[1];
rz(-1.0464484) q[1];
sx q[1];
rz(-1.4078543) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5124328) q[3];
sx q[3];
rz(-1.0608309) q[3];
sx q[3];
rz(2.8327018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87806988) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(0.44352201) q[2];
rz(-1.589132) q[3];
sx q[3];
rz(-0.3905206) q[3];
sx q[3];
rz(-2.922557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.012101128) q[0];
sx q[0];
rz(-0.87915593) q[0];
sx q[0];
rz(-2.7398859) q[0];
rz(-1.4178287) q[1];
sx q[1];
rz(-0.44770733) q[1];
sx q[1];
rz(1.1161425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6594769) q[0];
sx q[0];
rz(-1.2667313) q[0];
sx q[0];
rz(-0.67993865) q[0];
x q[1];
rz(-1.3858476) q[2];
sx q[2];
rz(-0.70025899) q[2];
sx q[2];
rz(1.8691693) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9520719) q[1];
sx q[1];
rz(-0.6743938) q[1];
sx q[1];
rz(2.2728069) q[1];
rz(-pi) q[2];
rz(-1.6468372) q[3];
sx q[3];
rz(-2.1180987) q[3];
sx q[3];
rz(0.28441498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22905722) q[2];
sx q[2];
rz(-2.8450862) q[2];
sx q[2];
rz(-2.1336446) q[2];
rz(-1.5063162) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(-0.87872163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85663831) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(2.304049) q[0];
rz(-1.8800927) q[1];
sx q[1];
rz(-1.3976169) q[1];
sx q[1];
rz(1.2417485) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1364131) q[0];
sx q[0];
rz(-1.8417999) q[0];
sx q[0];
rz(-2.6684929) q[0];
rz(-1.2707082) q[2];
sx q[2];
rz(-0.63020149) q[2];
sx q[2];
rz(0.84033191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82363331) q[1];
sx q[1];
rz(-1.9775864) q[1];
sx q[1];
rz(-0.49655668) q[1];
rz(-pi) q[2];
rz(0.99026525) q[3];
sx q[3];
rz(-1.8860236) q[3];
sx q[3];
rz(0.48993708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9813098) q[2];
sx q[2];
rz(-2.1872988) q[2];
sx q[2];
rz(-2.1295638) q[2];
rz(0.99803287) q[3];
sx q[3];
rz(-0.73603743) q[3];
sx q[3];
rz(1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0618458) q[0];
sx q[0];
rz(-0.80119067) q[0];
sx q[0];
rz(-0.27859846) q[0];
rz(-2.8246236) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(0.46151361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7636895) q[0];
sx q[0];
rz(-2.9660241) q[0];
sx q[0];
rz(-1.8493091) q[0];
rz(-0.95038173) q[2];
sx q[2];
rz(-1.912552) q[2];
sx q[2];
rz(-2.1555962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0089142) q[1];
sx q[1];
rz(-1.8092833) q[1];
sx q[1];
rz(-3.0452252) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0449681) q[3];
sx q[3];
rz(-1.0276252) q[3];
sx q[3];
rz(2.4718049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7757814) q[2];
sx q[2];
rz(-1.8041939) q[2];
sx q[2];
rz(0.27654761) q[2];
rz(2.5413399) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14715956) q[0];
sx q[0];
rz(-0.58013791) q[0];
sx q[0];
rz(-2.4603727) q[0];
rz(-0.45411202) q[1];
sx q[1];
rz(-1.9845767) q[1];
sx q[1];
rz(2.5195794) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2300782) q[0];
sx q[0];
rz(-2.921836) q[0];
sx q[0];
rz(0.4498555) q[0];
x q[1];
rz(2.3756681) q[2];
sx q[2];
rz(-2.3422547) q[2];
sx q[2];
rz(2.9098985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45251915) q[1];
sx q[1];
rz(-1.1645242) q[1];
sx q[1];
rz(-2.3157332) q[1];
rz(-pi) q[2];
rz(1.2219023) q[3];
sx q[3];
rz(-1.8233646) q[3];
sx q[3];
rz(-3.0608197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5155718) q[2];
sx q[2];
rz(-2.5280648) q[2];
sx q[2];
rz(-0.75508368) q[2];
rz(0.10609047) q[3];
sx q[3];
rz(-1.2664436) q[3];
sx q[3];
rz(2.6809926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945011) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(3.0754572) q[0];
rz(-3.0905837) q[1];
sx q[1];
rz(-0.74791932) q[1];
sx q[1];
rz(-1.8775108) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98718131) q[0];
sx q[0];
rz(-2.4307152) q[0];
sx q[0];
rz(2.1412503) q[0];
rz(-pi) q[1];
rz(-1.2602006) q[2];
sx q[2];
rz(-0.87970886) q[2];
sx q[2];
rz(0.6882918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6173827) q[1];
sx q[1];
rz(-2.2508989) q[1];
sx q[1];
rz(1.0587803) q[1];
rz(-0.73406666) q[3];
sx q[3];
rz(-2.5841004) q[3];
sx q[3];
rz(0.25171425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85822004) q[2];
sx q[2];
rz(-1.1274575) q[2];
sx q[2];
rz(-0.14632012) q[2];
rz(-0.60394168) q[3];
sx q[3];
rz(-0.54655176) q[3];
sx q[3];
rz(1.4124136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.716575) q[0];
sx q[0];
rz(-1.1973493) q[0];
sx q[0];
rz(-0.091751598) q[0];
rz(-0.07235202) q[1];
sx q[1];
rz(-1.3183343) q[1];
sx q[1];
rz(0.66973698) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2429072) q[0];
sx q[0];
rz(-1.4353485) q[0];
sx q[0];
rz(-0.29512914) q[0];
rz(-pi) q[1];
rz(1.4402499) q[2];
sx q[2];
rz(-1.1645082) q[2];
sx q[2];
rz(0.31995904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32423985) q[1];
sx q[1];
rz(-1.8756556) q[1];
sx q[1];
rz(-1.4190116) q[1];
x q[2];
rz(0.1454788) q[3];
sx q[3];
rz(-2.0091741) q[3];
sx q[3];
rz(-0.22989014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.79121315) q[2];
sx q[2];
rz(-2.4243441) q[2];
sx q[2];
rz(0.72163248) q[2];
rz(-2.3676938) q[3];
sx q[3];
rz(-2.1589203) q[3];
sx q[3];
rz(2.6017046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4197107) q[0];
sx q[0];
rz(-1.1996491) q[0];
sx q[0];
rz(2.6019959) q[0];
rz(2.3609912) q[1];
sx q[1];
rz(-1.8949948) q[1];
sx q[1];
rz(-0.4221198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2950384) q[0];
sx q[0];
rz(-2.1384412) q[0];
sx q[0];
rz(-0.66935434) q[0];
rz(0.11123734) q[2];
sx q[2];
rz(-0.49477494) q[2];
sx q[2];
rz(-0.24559427) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5752435) q[1];
sx q[1];
rz(-0.48108992) q[1];
sx q[1];
rz(1.4884557) q[1];
x q[2];
rz(-3.0738281) q[3];
sx q[3];
rz(-0.58971206) q[3];
sx q[3];
rz(0.36524352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.098794) q[2];
sx q[2];
rz(-2.1197987) q[2];
sx q[2];
rz(-2.5928024) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(-2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3720836) q[0];
sx q[0];
rz(-0.39396572) q[0];
sx q[0];
rz(-0.60419303) q[0];
rz(-1.6744042) q[1];
sx q[1];
rz(-1.3507651) q[1];
sx q[1];
rz(1.1473354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1579188) q[0];
sx q[0];
rz(-2.3566453) q[0];
sx q[0];
rz(-0.60529937) q[0];
rz(-pi) q[1];
rz(2.7962923) q[2];
sx q[2];
rz(-2.3757114) q[2];
sx q[2];
rz(1.6102546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5529788) q[1];
sx q[1];
rz(-2.0535894) q[1];
sx q[1];
rz(-1.5881513) q[1];
rz(1.9705521) q[3];
sx q[3];
rz(-1.2467524) q[3];
sx q[3];
rz(2.7017665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6818105) q[2];
sx q[2];
rz(-1.9777538) q[2];
sx q[2];
rz(-2.5414844) q[2];
rz(-0.72276226) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(2.7616937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4031274) q[0];
sx q[0];
rz(-1.6275788) q[0];
sx q[0];
rz(1.7049261) q[0];
rz(1.5557095) q[1];
sx q[1];
rz(-1.7441505) q[1];
sx q[1];
rz(2.2081262) q[1];
rz(0.97572648) q[2];
sx q[2];
rz(-2.3806934) q[2];
sx q[2];
rz(-1.0687547) q[2];
rz(0.71674552) q[3];
sx q[3];
rz(-2.8195753) q[3];
sx q[3];
rz(-1.4618518) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
