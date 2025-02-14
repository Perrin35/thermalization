OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1342993) q[0];
sx q[0];
rz(-1.1987885) q[0];
sx q[0];
rz(-0.78517908) q[0];
rz(2.5685318) q[1];
sx q[1];
rz(-0.95102024) q[1];
sx q[1];
rz(-2.5101488) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307237) q[0];
sx q[0];
rz(-2.0018303) q[0];
sx q[0];
rz(-2.2742989) q[0];
rz(-pi) q[1];
x q[1];
rz(2.417067) q[2];
sx q[2];
rz(-2.2093378) q[2];
sx q[2];
rz(-2.2341626) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50791477) q[1];
sx q[1];
rz(-0.79601015) q[1];
sx q[1];
rz(-0.62263864) q[1];
rz(2.9358125) q[3];
sx q[3];
rz(-0.45318402) q[3];
sx q[3];
rz(0.024621016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50420061) q[2];
sx q[2];
rz(-1.7039508) q[2];
sx q[2];
rz(-0.27153095) q[2];
rz(-0.046393752) q[3];
sx q[3];
rz(-1.4103187) q[3];
sx q[3];
rz(2.7739649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72535998) q[0];
sx q[0];
rz(-3.1144436) q[0];
sx q[0];
rz(-2.695988) q[0];
rz(-0.36000571) q[1];
sx q[1];
rz(-1.5574734) q[1];
sx q[1];
rz(-0.97703591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8603823) q[0];
sx q[0];
rz(-1.0458676) q[0];
sx q[0];
rz(2.1674744) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72190921) q[2];
sx q[2];
rz(-1.0470265) q[2];
sx q[2];
rz(-2.3350451) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9865446) q[1];
sx q[1];
rz(-0.98216559) q[1];
sx q[1];
rz(2.9595093) q[1];
rz(-pi) q[2];
rz(2.2350639) q[3];
sx q[3];
rz(-2.1593173) q[3];
sx q[3];
rz(1.9389467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52686024) q[2];
sx q[2];
rz(-1.7975668) q[2];
sx q[2];
rz(2.1993707) q[2];
rz(-0.49643907) q[3];
sx q[3];
rz(-1.4632635) q[3];
sx q[3];
rz(2.0739323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0838858) q[0];
sx q[0];
rz(-1.3297465) q[0];
sx q[0];
rz(2.4853793) q[0];
rz(0.48037848) q[1];
sx q[1];
rz(-0.96157938) q[1];
sx q[1];
rz(-0.61294714) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71343173) q[0];
sx q[0];
rz(-1.1640004) q[0];
sx q[0];
rz(-3.0359546) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0080994227) q[2];
sx q[2];
rz(-0.73609422) q[2];
sx q[2];
rz(-1.2978293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7483116) q[1];
sx q[1];
rz(-1.5419811) q[1];
sx q[1];
rz(-3.1080212) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5016236) q[3];
sx q[3];
rz(-2.0687447) q[3];
sx q[3];
rz(-0.70424265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2059242) q[2];
sx q[2];
rz(-1.0574477) q[2];
sx q[2];
rz(-2.617188) q[2];
rz(-2.2762401) q[3];
sx q[3];
rz(-1.0582346) q[3];
sx q[3];
rz(0.66044468) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48786369) q[0];
sx q[0];
rz(-1.8346584) q[0];
sx q[0];
rz(-3.0453239) q[0];
rz(0.013821566) q[1];
sx q[1];
rz(-2.6410069) q[1];
sx q[1];
rz(-1.0523419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88599624) q[0];
sx q[0];
rz(-1.692472) q[0];
sx q[0];
rz(-0.047329655) q[0];
rz(-pi) q[1];
rz(2.0180447) q[2];
sx q[2];
rz(-2.5442225) q[2];
sx q[2];
rz(1.1505466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87409243) q[1];
sx q[1];
rz(-2.7639031) q[1];
sx q[1];
rz(-0.51376359) q[1];
x q[2];
rz(-2.0963444) q[3];
sx q[3];
rz(-1.6017672) q[3];
sx q[3];
rz(-0.82079673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4200165) q[2];
sx q[2];
rz(-0.53109157) q[2];
sx q[2];
rz(-2.6225923) q[2];
rz(-2.0429933) q[3];
sx q[3];
rz(-1.4562621) q[3];
sx q[3];
rz(2.4115653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5191583) q[0];
sx q[0];
rz(-2.9485478) q[0];
sx q[0];
rz(0.73424196) q[0];
rz(1.2892067) q[1];
sx q[1];
rz(-1.0451008) q[1];
sx q[1];
rz(-2.2330914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8216128) q[0];
sx q[0];
rz(-0.93324086) q[0];
sx q[0];
rz(2.430401) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5703326) q[2];
sx q[2];
rz(-1.6831584) q[2];
sx q[2];
rz(-1.4552356) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1724194) q[1];
sx q[1];
rz(-2.4052019) q[1];
sx q[1];
rz(-1.4672155) q[1];
rz(1.8374422) q[3];
sx q[3];
rz(-1.5943267) q[3];
sx q[3];
rz(-1.8950099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.044540731) q[2];
sx q[2];
rz(-2.3553607) q[2];
sx q[2];
rz(2.457974) q[2];
rz(-2.0475856) q[3];
sx q[3];
rz(-1.8180327) q[3];
sx q[3];
rz(-1.5692284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.053726824) q[0];
sx q[0];
rz(-1.6931067) q[0];
sx q[0];
rz(2.7159765) q[0];
rz(0.60676891) q[1];
sx q[1];
rz(-0.68777045) q[1];
sx q[1];
rz(-1.9389796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9894597) q[0];
sx q[0];
rz(-2.0589863) q[0];
sx q[0];
rz(-0.91849416) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21738217) q[2];
sx q[2];
rz(-0.067408407) q[2];
sx q[2];
rz(1.6081152) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6708199) q[1];
sx q[1];
rz(-2.2774376) q[1];
sx q[1];
rz(0.67292893) q[1];
rz(-pi) q[2];
rz(0.7814066) q[3];
sx q[3];
rz(-2.8202116) q[3];
sx q[3];
rz(1.7226294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46835029) q[2];
sx q[2];
rz(-0.52933669) q[2];
sx q[2];
rz(0.62874111) q[2];
rz(2.8041503) q[3];
sx q[3];
rz(-1.3597666) q[3];
sx q[3];
rz(0.76702816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5479946) q[0];
sx q[0];
rz(-0.070040528) q[0];
sx q[0];
rz(-1.3399757) q[0];
rz(-0.10546671) q[1];
sx q[1];
rz(-2.1601845) q[1];
sx q[1];
rz(-2.9973105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7026065) q[0];
sx q[0];
rz(-2.6487338) q[0];
sx q[0];
rz(-0.65123691) q[0];
rz(-pi) q[1];
rz(0.86267306) q[2];
sx q[2];
rz(-1.4221592) q[2];
sx q[2];
rz(1.7903863) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0622341) q[1];
sx q[1];
rz(-2.7309504) q[1];
sx q[1];
rz(-2.5985318) q[1];
x q[2];
rz(0.67146164) q[3];
sx q[3];
rz(-1.1290159) q[3];
sx q[3];
rz(-2.9908671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1062801) q[2];
sx q[2];
rz(-1.148618) q[2];
sx q[2];
rz(2.3036352) q[2];
rz(2.641814) q[3];
sx q[3];
rz(-0.79334799) q[3];
sx q[3];
rz(-0.70118538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226294) q[0];
sx q[0];
rz(-2.0327649) q[0];
sx q[0];
rz(0.75123373) q[0];
rz(-1.2535837) q[1];
sx q[1];
rz(-1.3026078) q[1];
sx q[1];
rz(-1.368103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.591044) q[0];
sx q[0];
rz(-0.97752377) q[0];
sx q[0];
rz(0.53925487) q[0];
x q[1];
rz(1.8644731) q[2];
sx q[2];
rz(-1.2674567) q[2];
sx q[2];
rz(-1.9544698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3093395) q[1];
sx q[1];
rz(-2.3414438) q[1];
sx q[1];
rz(0.11804987) q[1];
x q[2];
rz(-0.99686025) q[3];
sx q[3];
rz(-1.878345) q[3];
sx q[3];
rz(2.779863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83296835) q[2];
sx q[2];
rz(-0.37165752) q[2];
sx q[2];
rz(0.0011681636) q[2];
rz(2.809281) q[3];
sx q[3];
rz(-1.6806335) q[3];
sx q[3];
rz(2.6619679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7697656) q[0];
sx q[0];
rz(-1.2170987) q[0];
sx q[0];
rz(1.0960854) q[0];
rz(-1.7425884) q[1];
sx q[1];
rz(-1.5053446) q[1];
sx q[1];
rz(-0.91986626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8815589) q[0];
sx q[0];
rz(-1.1333916) q[0];
sx q[0];
rz(-3.0946671) q[0];
rz(-pi) q[1];
rz(1.9894753) q[2];
sx q[2];
rz(-1.8279462) q[2];
sx q[2];
rz(-0.60794431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6429872) q[1];
sx q[1];
rz(-0.66434331) q[1];
sx q[1];
rz(-0.24533116) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9971531) q[3];
sx q[3];
rz(-1.0424992) q[3];
sx q[3];
rz(-2.4424255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9225191) q[2];
sx q[2];
rz(-0.24524958) q[2];
sx q[2];
rz(-1.6314487) q[2];
rz(2.505693) q[3];
sx q[3];
rz(-2.4363775) q[3];
sx q[3];
rz(-2.6270134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9566327) q[0];
sx q[0];
rz(-2.2745467) q[0];
sx q[0];
rz(1.6459203) q[0];
rz(3.0236859) q[1];
sx q[1];
rz(-1.3943744) q[1];
sx q[1];
rz(-0.25679055) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816329) q[0];
sx q[0];
rz(-1.5149533) q[0];
sx q[0];
rz(1.1373489) q[0];
rz(-pi) q[1];
rz(-0.43440993) q[2];
sx q[2];
rz(-2.5108178) q[2];
sx q[2];
rz(1.3303016) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3658095) q[1];
sx q[1];
rz(-1.3089615) q[1];
sx q[1];
rz(-0.29579137) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.58415) q[3];
sx q[3];
rz(-0.62288219) q[3];
sx q[3];
rz(-2.0940144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7439277) q[2];
sx q[2];
rz(-2.3981018) q[2];
sx q[2];
rz(2.1916154) q[2];
rz(2.5395565) q[3];
sx q[3];
rz(-1.2369912) q[3];
sx q[3];
rz(-0.33815798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.4468741) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(1.0097722) q[1];
sx q[1];
rz(-0.92242454) q[1];
sx q[1];
rz(-0.19933272) q[1];
rz(1.2620325) q[2];
sx q[2];
rz(-1.2876576) q[2];
sx q[2];
rz(1.8907945) q[2];
rz(-1.1789523) q[3];
sx q[3];
rz(-2.3096991) q[3];
sx q[3];
rz(0.25326412) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
