OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1256975) q[0];
sx q[0];
rz(-0.96244922) q[0];
sx q[0];
rz(-0.35981052) q[0];
x q[1];
rz(1.4719226) q[2];
sx q[2];
rz(-0.31497248) q[2];
sx q[2];
rz(-2.0801983) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4830556) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(0.87941054) q[1];
x q[2];
rz(2.764774) q[3];
sx q[3];
rz(-1.0816649) q[3];
sx q[3];
rz(0.34987846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(-0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(-0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(-0.91631779) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57740649) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(2.2400411) q[0];
rz(-pi) q[1];
rz(1.4257405) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(-1.8650101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9333222) q[1];
sx q[1];
rz(-1.2332488) q[1];
sx q[1];
rz(2.7949105) q[1];
x q[2];
rz(-1.5171492) q[3];
sx q[3];
rz(-0.39857736) q[3];
sx q[3];
rz(-0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(2.6170513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0073111) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(1.4071464) q[0];
rz(-pi) q[1];
rz(-1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(1.7922572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14568612) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(-1.5476336) q[1];
rz(-pi) q[2];
rz(2.9933661) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(1.9849329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77164578) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(0.088949732) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3084429) q[0];
sx q[0];
rz(-1.2763378) q[0];
sx q[0];
rz(-1.4008646) q[0];
rz(-pi) q[1];
rz(-1.3454516) q[2];
sx q[2];
rz(-2.1225404) q[2];
sx q[2];
rz(1.6793041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.888962) q[1];
sx q[1];
rz(-1.8969715) q[1];
sx q[1];
rz(2.3701282) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(-1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(2.518667) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.6437644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9153584) q[0];
sx q[0];
rz(-2.3225473) q[0];
sx q[0];
rz(-2.782343) q[0];
x q[1];
rz(3.0450902) q[2];
sx q[2];
rz(-2.3415903) q[2];
sx q[2];
rz(2.2068791) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0651107) q[1];
sx q[1];
rz(-1.5434192) q[1];
sx q[1];
rz(-0.70365023) q[1];
rz(-pi) q[2];
rz(0.51550128) q[3];
sx q[3];
rz(-1.4294129) q[3];
sx q[3];
rz(-0.055671234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(1.4354338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88046056) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(2.3418505) q[0];
rz(-pi) q[1];
rz(2.1204505) q[2];
sx q[2];
rz(-1.5205985) q[2];
sx q[2];
rz(-0.93027885) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3458851) q[1];
sx q[1];
rz(-0.42783034) q[1];
sx q[1];
rz(1.8756299) q[1];
x q[2];
rz(1.3330323) q[3];
sx q[3];
rz(-1.5441582) q[3];
sx q[3];
rz(-2.3029072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(1.2795992) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(2.1320027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44430915) q[0];
sx q[0];
rz(-0.59519207) q[0];
sx q[0];
rz(-2.5213581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7570653) q[2];
sx q[2];
rz(-1.782801) q[2];
sx q[2];
rz(2.7146313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.13751444) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(-0.60933463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3791802) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(1.4095969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(1.002731) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49884477) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(-2.902466) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(-0.4447287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0895773) q[0];
sx q[0];
rz(-1.9964295) q[0];
sx q[0];
rz(2.1968967) q[0];
rz(-pi) q[1];
rz(0.93034805) q[2];
sx q[2];
rz(-1.010251) q[2];
sx q[2];
rz(2.408574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4152894) q[1];
sx q[1];
rz(-0.88978926) q[1];
sx q[1];
rz(2.9773832) q[1];
rz(-pi) q[2];
rz(-2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5765604) q[3];
sx q[3];
rz(-1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(0.28717336) q[0];
rz(-0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(0.33219355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9438292) q[0];
sx q[0];
rz(-1.8837067) q[0];
sx q[0];
rz(2.3373332) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(-0.38682129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6576782) q[1];
sx q[1];
rz(-0.49424833) q[1];
sx q[1];
rz(-2.5792522) q[1];
rz(-pi) q[2];
rz(2.6094749) q[3];
sx q[3];
rz(-1.2174774) q[3];
sx q[3];
rz(1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-0.65657842) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(-0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018325018) q[0];
sx q[0];
rz(-1.3402481) q[0];
sx q[0];
rz(0.075320764) q[0];
x q[1];
rz(2.6539831) q[2];
sx q[2];
rz(-1.3541823) q[2];
sx q[2];
rz(-0.85607869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0925797) q[1];
sx q[1];
rz(-1.6515459) q[1];
sx q[1];
rz(0.32596522) q[1];
rz(-1.6700527) q[3];
sx q[3];
rz(-2.0072862) q[3];
sx q[3];
rz(1.0222767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142674) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-1.1709471) q[2];
sx q[2];
rz(-0.56213899) q[2];
sx q[2];
rz(3.0664372) q[2];
rz(-2.1847235) q[3];
sx q[3];
rz(-1.3812243) q[3];
sx q[3];
rz(0.45637043) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
