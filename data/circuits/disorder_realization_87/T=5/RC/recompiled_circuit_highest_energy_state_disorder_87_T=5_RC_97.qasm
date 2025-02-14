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
rz(-1.0083899) q[1];
sx q[1];
rz(-2.402306) q[1];
sx q[1];
rz(-0.33831212) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9341921) q[0];
sx q[0];
rz(-1.3284042) q[0];
sx q[0];
rz(1.5336179) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9954483) q[2];
sx q[2];
rz(-2.3156347) q[2];
sx q[2];
rz(0.34817696) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4614726) q[1];
sx q[1];
rz(-1.0306038) q[1];
sx q[1];
rz(-1.4137181) q[1];
rz(-1.385072) q[3];
sx q[3];
rz(-0.76558569) q[3];
sx q[3];
rz(3.035745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7484635) q[2];
sx q[2];
rz(-1.7270361) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2748134) q[0];
sx q[0];
rz(-1.4115189) q[0];
sx q[0];
rz(-2.859512) q[0];
rz(-0.2433978) q[1];
sx q[1];
rz(-1.2671821) q[1];
sx q[1];
rz(2.3928221) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7215749) q[0];
sx q[0];
rz(-0.87792464) q[0];
sx q[0];
rz(1.3189032) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1339945) q[2];
sx q[2];
rz(-1.6924592) q[2];
sx q[2];
rz(1.0341782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91606027) q[1];
sx q[1];
rz(-0.5114564) q[1];
sx q[1];
rz(1.3813853) q[1];
x q[2];
rz(-1.5010819) q[3];
sx q[3];
rz(-1.3718714) q[3];
sx q[3];
rz(0.45123842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8613209) q[2];
sx q[2];
rz(-1.6161641) q[2];
sx q[2];
rz(-1.3410428) q[2];
rz(-1.1567814) q[3];
sx q[3];
rz(-2.3564434) q[3];
sx q[3];
rz(2.3780499) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3007616) q[0];
sx q[0];
rz(-0.16591993) q[0];
sx q[0];
rz(-2.4842343) q[0];
rz(1.1454469) q[1];
sx q[1];
rz(-0.95445389) q[1];
sx q[1];
rz(2.3232536) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377146) q[0];
sx q[0];
rz(-2.7174207) q[0];
sx q[0];
rz(0.61221497) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4253699) q[2];
sx q[2];
rz(-2.8231695) q[2];
sx q[2];
rz(2.2216036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.9200478) q[1];
sx q[1];
rz(-0.86600477) q[1];
sx q[1];
rz(-0.66001604) q[1];
rz(-3.0679296) q[3];
sx q[3];
rz(-2.0233002) q[3];
sx q[3];
rz(-1.0301925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98693097) q[2];
sx q[2];
rz(-1.8795452) q[2];
sx q[2];
rz(1.3531125) q[2];
rz(2.1766369) q[3];
sx q[3];
rz(-1.8392287) q[3];
sx q[3];
rz(-0.42199782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593909) q[0];
sx q[0];
rz(-2.4762479) q[0];
sx q[0];
rz(-1.9406142) q[0];
rz(3.1383842) q[1];
sx q[1];
rz(-1.708958) q[1];
sx q[1];
rz(-2.9046955) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2819716) q[0];
sx q[0];
rz(-1.0399204) q[0];
sx q[0];
rz(-2.6174699) q[0];
rz(0.63726823) q[2];
sx q[2];
rz(-1.6760525) q[2];
sx q[2];
rz(0.95544514) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.264852) q[1];
sx q[1];
rz(-2.5036466) q[1];
sx q[1];
rz(1.2400342) q[1];
rz(-pi) q[2];
rz(-2.3637216) q[3];
sx q[3];
rz(-0.98613769) q[3];
sx q[3];
rz(-2.8608222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0995471) q[2];
sx q[2];
rz(-1.8501661) q[2];
sx q[2];
rz(0.42638865) q[2];
rz(2.1060627) q[3];
sx q[3];
rz(-1.6798881) q[3];
sx q[3];
rz(0.6216014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34430382) q[0];
sx q[0];
rz(-3.099589) q[0];
sx q[0];
rz(2.7916743) q[0];
rz(-1.2087076) q[1];
sx q[1];
rz(-1.8588926) q[1];
sx q[1];
rz(0.0016317687) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7662273) q[0];
sx q[0];
rz(-1.3266801) q[0];
sx q[0];
rz(1.7422416) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6095743) q[2];
sx q[2];
rz(-1.8680405) q[2];
sx q[2];
rz(-2.4448038) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35910159) q[1];
sx q[1];
rz(-1.3939855) q[1];
sx q[1];
rz(-0.048080877) q[1];
x q[2];
rz(-2.2285813) q[3];
sx q[3];
rz(-0.33644852) q[3];
sx q[3];
rz(-2.2483765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95167595) q[2];
sx q[2];
rz(-2.2168171) q[2];
sx q[2];
rz(2.9803989) q[2];
rz(-2.7023884) q[3];
sx q[3];
rz(-0.74857155) q[3];
sx q[3];
rz(-2.2883033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31140232) q[0];
sx q[0];
rz(-1.4051733) q[0];
sx q[0];
rz(0.79469529) q[0];
rz(1.7757802) q[1];
sx q[1];
rz(-2.3410485) q[1];
sx q[1];
rz(-0.47439233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2524388) q[0];
sx q[0];
rz(-1.948101) q[0];
sx q[0];
rz(-1.0214367) q[0];
rz(-pi) q[1];
rz(2.5849336) q[2];
sx q[2];
rz(-1.8881329) q[2];
sx q[2];
rz(-0.59653004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77955708) q[1];
sx q[1];
rz(-1.0048702) q[1];
sx q[1];
rz(2.0574244) q[1];
x q[2];
rz(2.1204409) q[3];
sx q[3];
rz(-2.4826038) q[3];
sx q[3];
rz(-0.07168183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5279493) q[2];
sx q[2];
rz(-1.2443292) q[2];
sx q[2];
rz(-2.6141686) q[2];
rz(0.6066277) q[3];
sx q[3];
rz(-0.18692034) q[3];
sx q[3];
rz(-2.0691779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793176) q[0];
sx q[0];
rz(-0.19392218) q[0];
sx q[0];
rz(-0.79757565) q[0];
rz(2.2480615) q[1];
sx q[1];
rz(-2.306566) q[1];
sx q[1];
rz(0.28018793) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55927568) q[0];
sx q[0];
rz(-1.7452551) q[0];
sx q[0];
rz(-1.1691537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75887078) q[2];
sx q[2];
rz(-0.810383) q[2];
sx q[2];
rz(-0.23960613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5366027) q[1];
sx q[1];
rz(-1.0883253) q[1];
sx q[1];
rz(-3.1274113) q[1];
x q[2];
rz(-1.0605162) q[3];
sx q[3];
rz(-2.4852537) q[3];
sx q[3];
rz(2.9440232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4013275) q[2];
sx q[2];
rz(-1.836931) q[2];
sx q[2];
rz(-1.9291482) q[2];
rz(-0.23981833) q[3];
sx q[3];
rz(-2.4739154) q[3];
sx q[3];
rz(-0.56813204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577393) q[0];
sx q[0];
rz(-1.4171866) q[0];
sx q[0];
rz(-1.3215815) q[0];
rz(2.9603738) q[1];
sx q[1];
rz(-1.8099338) q[1];
sx q[1];
rz(3.0785353) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4806594) q[0];
sx q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(0.78137915) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8562154) q[2];
sx q[2];
rz(-1.4865424) q[2];
sx q[2];
rz(-2.7610996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5279505) q[1];
sx q[1];
rz(-2.7220352) q[1];
sx q[1];
rz(-2.502524) q[1];
x q[2];
rz(-1.5531055) q[3];
sx q[3];
rz(-2.9790386) q[3];
sx q[3];
rz(1.7944195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41398373) q[2];
sx q[2];
rz(-2.0435645) q[2];
sx q[2];
rz(-1.6723527) q[2];
rz(-1.3937048) q[3];
sx q[3];
rz(-1.4398451) q[3];
sx q[3];
rz(0.20974717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.1446153) q[0];
rz(0.793055) q[1];
sx q[1];
rz(-1.4184364) q[1];
sx q[1];
rz(1.1096032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5376268) q[0];
sx q[0];
rz(-1.0111799) q[0];
sx q[0];
rz(1.9163314) q[0];
rz(-pi) q[1];
rz(-2.6656261) q[2];
sx q[2];
rz(-2.360811) q[2];
sx q[2];
rz(2.4300574) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3396757) q[1];
sx q[1];
rz(-1.1028506) q[1];
sx q[1];
rz(-2.1956594) q[1];
x q[2];
rz(0.22623541) q[3];
sx q[3];
rz(-2.492815) q[3];
sx q[3];
rz(-1.7895376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7290466) q[2];
sx q[2];
rz(-1.6763326) q[2];
sx q[2];
rz(-1.868978) q[2];
rz(-2.0160969) q[3];
sx q[3];
rz(-0.19386217) q[3];
sx q[3];
rz(-0.079843609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.37128714) q[0];
sx q[0];
rz(-1.9885539) q[0];
sx q[0];
rz(0.09819296) q[0];
rz(1.498361) q[1];
sx q[1];
rz(-1.1474835) q[1];
sx q[1];
rz(0.8383382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22238092) q[0];
sx q[0];
rz(-1.7475238) q[0];
sx q[0];
rz(-2.3803898) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8894601) q[2];
sx q[2];
rz(-1.6246535) q[2];
sx q[2];
rz(-2.5359652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34217087) q[1];
sx q[1];
rz(-1.2890745) q[1];
sx q[1];
rz(-1.1200302) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7097852) q[3];
sx q[3];
rz(-1.7718959) q[3];
sx q[3];
rz(-2.3922298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3457634) q[2];
sx q[2];
rz(-2.2553307) q[2];
sx q[2];
rz(-0.32540992) q[2];
rz(2.365153) q[3];
sx q[3];
rz(-1.0151981) q[3];
sx q[3];
rz(-2.5344892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0008739) q[0];
sx q[0];
rz(-1.9372531) q[0];
sx q[0];
rz(-1.0304864) q[0];
rz(-0.28618947) q[1];
sx q[1];
rz(-1.9356526) q[1];
sx q[1];
rz(-2.0331358) q[1];
rz(1.1775511) q[2];
sx q[2];
rz(-2.2920592) q[2];
sx q[2];
rz(-0.68790676) q[2];
rz(-1.9097044) q[3];
sx q[3];
rz(-1.00885) q[3];
sx q[3];
rz(2.5918617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
