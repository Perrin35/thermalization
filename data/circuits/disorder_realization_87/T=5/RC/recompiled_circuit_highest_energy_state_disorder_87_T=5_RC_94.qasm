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
rz(5.7558007) q[0];
sx q[0];
rz(5.7100073) q[0];
sx q[0];
rz(10.041458) q[0];
rz(2.1332027) q[1];
sx q[1];
rz(-0.73928666) q[1];
sx q[1];
rz(0.33831212) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341921) q[0];
sx q[0];
rz(-1.8131885) q[0];
sx q[0];
rz(1.5336179) q[0];
rz(-pi) q[1];
rz(1.1461444) q[2];
sx q[2];
rz(-2.3156347) q[2];
sx q[2];
rz(-2.7934157) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4614726) q[1];
sx q[1];
rz(-1.0306038) q[1];
sx q[1];
rz(1.7278746) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1756535) q[3];
sx q[3];
rz(-2.3199816) q[3];
sx q[3];
rz(2.780811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7484635) q[2];
sx q[2];
rz(-1.4145565) q[2];
sx q[2];
rz(2.4836922) q[2];
rz(0.45514485) q[3];
sx q[3];
rz(-0.15076605) q[3];
sx q[3];
rz(-1.3078825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748134) q[0];
sx q[0];
rz(-1.7300737) q[0];
sx q[0];
rz(-0.28208062) q[0];
rz(2.8981949) q[1];
sx q[1];
rz(-1.8744105) q[1];
sx q[1];
rz(0.74877053) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42001777) q[0];
sx q[0];
rz(-0.87792464) q[0];
sx q[0];
rz(1.3189032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.44913) q[2];
sx q[2];
rz(-1.5632544) q[2];
sx q[2];
rz(2.6040524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1324596) q[1];
sx q[1];
rz(-1.0693502) q[1];
sx q[1];
rz(-0.10528586) q[1];
x q[2];
rz(-0.19939662) q[3];
sx q[3];
rz(-1.5024589) q[3];
sx q[3];
rz(2.0082366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28027174) q[2];
sx q[2];
rz(-1.5254285) q[2];
sx q[2];
rz(-1.8005499) q[2];
rz(1.9848112) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(-0.7635428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84083104) q[0];
sx q[0];
rz(-0.16591993) q[0];
sx q[0];
rz(2.4842343) q[0];
rz(1.9961458) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(0.81833902) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3439804) q[0];
sx q[0];
rz(-1.8095785) q[0];
sx q[0];
rz(0.35399951) q[0];
rz(-pi) q[1];
rz(-3.0938593) q[2];
sx q[2];
rz(-1.2558508) q[2];
sx q[2];
rz(-0.76698179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.044825252) q[1];
sx q[1];
rz(-0.92508537) q[1];
sx q[1];
rz(2.195408) q[1];
rz(0.07366304) q[3];
sx q[3];
rz(-2.0233002) q[3];
sx q[3];
rz(2.1114001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1546617) q[2];
sx q[2];
rz(-1.2620474) q[2];
sx q[2];
rz(1.7884802) q[2];
rz(-0.96495572) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(2.7195948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956536) q[0];
sx q[0];
rz(-0.66534477) q[0];
sx q[0];
rz(1.9406142) q[0];
rz(3.1383842) q[1];
sx q[1];
rz(-1.708958) q[1];
sx q[1];
rz(-2.9046955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1457082) q[0];
sx q[0];
rz(-2.0170324) q[0];
sx q[0];
rz(-2.1666885) q[0];
x q[1];
rz(1.4400993) q[2];
sx q[2];
rz(-0.937619) q[2];
sx q[2];
rz(-2.6038632) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.264852) q[1];
sx q[1];
rz(-2.5036466) q[1];
sx q[1];
rz(1.2400342) q[1];
rz(-pi) q[2];
rz(-0.77787106) q[3];
sx q[3];
rz(-0.98613769) q[3];
sx q[3];
rz(-0.28077048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0995471) q[2];
sx q[2];
rz(-1.8501661) q[2];
sx q[2];
rz(0.42638865) q[2];
rz(1.03553) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(2.5199913) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34430382) q[0];
sx q[0];
rz(-0.042003691) q[0];
sx q[0];
rz(-0.34991831) q[0];
rz(1.2087076) q[1];
sx q[1];
rz(-1.2827001) q[1];
sx q[1];
rz(0.0016317687) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24626075) q[0];
sx q[0];
rz(-0.29731942) q[0];
sx q[0];
rz(2.5410482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29745488) q[2];
sx q[2];
rz(-1.5337197) q[2];
sx q[2];
rz(2.2789479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35910159) q[1];
sx q[1];
rz(-1.7476071) q[1];
sx q[1];
rz(-3.0935118) q[1];
rz(-pi) q[2];
rz(-1.3007845) q[3];
sx q[3];
rz(-1.367566) q[3];
sx q[3];
rz(-3.0940987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95167595) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(2.9803989) q[2];
rz(2.7023884) q[3];
sx q[3];
rz(-2.3930211) q[3];
sx q[3];
rz(0.85328931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31140232) q[0];
sx q[0];
rz(-1.4051733) q[0];
sx q[0];
rz(2.3468974) q[0];
rz(-1.3658124) q[1];
sx q[1];
rz(-0.8005442) q[1];
sx q[1];
rz(-2.6672003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5401986) q[0];
sx q[0];
rz(-2.0776333) q[0];
sx q[0];
rz(2.7066134) q[0];
rz(-pi) q[1];
rz(-2.5854255) q[2];
sx q[2];
rz(-2.5092297) q[2];
sx q[2];
rz(-1.4390505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51481828) q[1];
sx q[1];
rz(-1.1650135) q[1];
sx q[1];
rz(2.5184513) q[1];
x q[2];
rz(0.38444774) q[3];
sx q[3];
rz(-1.0214503) q[3];
sx q[3];
rz(0.58754301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61364335) q[2];
sx q[2];
rz(-1.2443292) q[2];
sx q[2];
rz(2.6141686) q[2];
rz(-2.534965) q[3];
sx q[3];
rz(-2.9546723) q[3];
sx q[3];
rz(2.0691779) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793176) q[0];
sx q[0];
rz(-2.9476705) q[0];
sx q[0];
rz(-2.344017) q[0];
rz(-2.2480615) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(2.8614047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2036616) q[0];
sx q[0];
rz(-1.1755921) q[0];
sx q[0];
rz(0.18919887) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1970388) q[2];
sx q[2];
rz(-2.124386) q[2];
sx q[2];
rz(1.9596726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.97238648) q[1];
sx q[1];
rz(-1.5833588) q[1];
sx q[1];
rz(-2.0533086) q[1];
rz(-pi) q[2];
rz(-2.0810764) q[3];
sx q[3];
rz(-2.4852537) q[3];
sx q[3];
rz(-2.9440232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4013275) q[2];
sx q[2];
rz(-1.3046616) q[2];
sx q[2];
rz(1.9291482) q[2];
rz(-2.9017743) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(-2.5734606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3577393) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(1.3215815) q[0];
rz(2.9603738) q[1];
sx q[1];
rz(-1.3316589) q[1];
sx q[1];
rz(-3.0785353) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94751271) q[0];
sx q[0];
rz(-0.79160684) q[0];
sx q[0];
rz(-2.9400154) q[0];
rz(-3.0538043) q[2];
sx q[2];
rz(-1.8551747) q[2];
sx q[2];
rz(1.926601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5279505) q[1];
sx q[1];
rz(-2.7220352) q[1];
sx q[1];
rz(-2.502524) q[1];
x q[2];
rz(-1.4082673) q[3];
sx q[3];
rz(-1.5679334) q[3];
sx q[3];
rz(-0.2410808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(1.6723527) q[2];
rz(-1.7478878) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(-0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53836981) q[0];
sx q[0];
rz(-2.5312238) q[0];
sx q[0];
rz(0.99697733) q[0];
rz(0.793055) q[1];
sx q[1];
rz(-1.7231562) q[1];
sx q[1];
rz(2.0319895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1556517) q[0];
sx q[0];
rz(-1.8619259) q[0];
sx q[0];
rz(2.554214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47596653) q[2];
sx q[2];
rz(-2.360811) q[2];
sx q[2];
rz(0.7115353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0834137) q[1];
sx q[1];
rz(-1.021487) q[1];
sx q[1];
rz(-0.55725248) q[1];
rz(-pi) q[2];
rz(-1.7392735) q[3];
sx q[3];
rz(-2.2003897) q[3];
sx q[3];
rz(1.633267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41254607) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(1.2726146) q[2];
rz(-2.0160969) q[3];
sx q[3];
rz(-2.9477305) q[3];
sx q[3];
rz(0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7703055) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(0.09819296) q[0];
rz(1.498361) q[1];
sx q[1];
rz(-1.9941092) q[1];
sx q[1];
rz(-0.8383382) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22238092) q[0];
sx q[0];
rz(-1.3940689) q[0];
sx q[0];
rz(0.76120283) q[0];
x q[1];
rz(-0.25213253) q[2];
sx q[2];
rz(-1.6246535) q[2];
sx q[2];
rz(0.60562741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0467207) q[1];
sx q[1];
rz(-1.1390242) q[1];
sx q[1];
rz(-2.8304965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3499979) q[3];
sx q[3];
rz(-1.9933369) q[3];
sx q[3];
rz(2.411946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7958293) q[2];
sx q[2];
rz(-0.886262) q[2];
sx q[2];
rz(-2.8161827) q[2];
rz(2.365153) q[3];
sx q[3];
rz(-2.1263945) q[3];
sx q[3];
rz(-0.6071035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0008739) q[0];
sx q[0];
rz(-1.2043395) q[0];
sx q[0];
rz(2.1111063) q[0];
rz(-2.8554032) q[1];
sx q[1];
rz(-1.20594) q[1];
sx q[1];
rz(1.1084569) q[1];
rz(0.4109702) q[2];
sx q[2];
rz(-2.3373418) q[2];
sx q[2];
rz(1.892754) q[2];
rz(0.48578942) q[3];
sx q[3];
rz(-0.64668568) q[3];
sx q[3];
rz(0.034737094) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
