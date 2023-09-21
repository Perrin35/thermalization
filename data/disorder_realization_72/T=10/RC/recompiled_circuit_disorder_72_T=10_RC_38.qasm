OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141657) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(-0.27557296) q[0];
rz(0.17272858) q[2];
sx q[2];
rz(-0.85731259) q[2];
sx q[2];
rz(-1.8634862) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.7826739) q[1];
sx q[1];
rz(-2.3896857) q[1];
sx q[1];
rz(-0.2207252) q[1];
rz(-pi) q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-3.0348274) q[3];
sx q[3];
rz(-1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-2.1710158) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(-2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-2.1388334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402055) q[0];
sx q[0];
rz(-1.1367072) q[0];
sx q[0];
rz(2.4853404) q[0];
rz(1.7582714) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(2.1263009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7748389) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(1.2971836) q[1];
x q[2];
rz(2.6456804) q[3];
sx q[3];
rz(-0.65973385) q[3];
sx q[3];
rz(2.0663313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(2.2580106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6813864) q[0];
sx q[0];
rz(-1.1707414) q[0];
sx q[0];
rz(-1.3921188) q[0];
rz(-1.9687037) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(0.18323252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2391501) q[1];
sx q[1];
rz(-2.855774) q[1];
sx q[1];
rz(1.6509389) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9900842) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-0.5501737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4470091) q[0];
sx q[0];
rz(-2.3059079) q[0];
sx q[0];
rz(-1.6741333) q[0];
rz(-pi) q[1];
rz(0.29454622) q[2];
sx q[2];
rz(-1.4809594) q[2];
sx q[2];
rz(1.3603269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(-2.1469841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6226193) q[3];
sx q[3];
rz(-1.7997777) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(-0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-0.91526389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(1.4907452) q[0];
x q[1];
rz(0.5251685) q[2];
sx q[2];
rz(-1.6095973) q[2];
sx q[2];
rz(-2.2280488) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3770959) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(1.7736969) q[1];
rz(2.3239273) q[3];
sx q[3];
rz(-0.92901232) q[3];
sx q[3];
rz(1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-2.6780224) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298033) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(-2.0407709) q[0];
rz(-pi) q[1];
rz(2.6655212) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(-1.6518041) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0551758) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(2.7820884) q[3];
sx q[3];
rz(-2.5264969) q[3];
sx q[3];
rz(-2.2744327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5863824) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(0.82715183) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92857498) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(-1.6638882) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0625213) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(-1.6342271) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9019482) q[1];
sx q[1];
rz(-2.6522954) q[1];
sx q[1];
rz(2.7326665) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(-1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(2.1408634) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(0.032756068) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(-2.2420922) q[0];
rz(2.6997386) q[2];
sx q[2];
rz(-1.0519069) q[2];
sx q[2];
rz(-1.9372802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36383648) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(1.9883518) q[1];
rz(-3.0524391) q[3];
sx q[3];
rz(-1.1925863) q[3];
sx q[3];
rz(-2.0034727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(2.5194871) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(-0.73293066) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-2.656235) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5068685) q[0];
sx q[0];
rz(-1.6244495) q[0];
sx q[0];
rz(3.0070478) q[0];
x q[1];
rz(0.23156667) q[2];
sx q[2];
rz(-0.61083691) q[2];
sx q[2];
rz(-1.6419186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(-2.0054714) q[1];
rz(-pi) q[2];
rz(-2.0210578) q[3];
sx q[3];
rz(-2.1053208) q[3];
sx q[3];
rz(0.083732097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(1.7074701) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(0.5982582) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.7822595) q[0];
rz(-pi) q[1];
rz(1.7190476) q[2];
sx q[2];
rz(-1.3607927) q[2];
sx q[2];
rz(-1.6809747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1640537) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(-1.8305199) q[1];
rz(-pi) q[2];
rz(1.7353021) q[3];
sx q[3];
rz(-1.4146155) q[3];
sx q[3];
rz(-0.7848878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(2.4895978) q[2];
rz(-2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-0.61959735) q[2];
sx q[2];
rz(-1.4705114) q[2];
sx q[2];
rz(0.32497139) q[2];
rz(0.89647722) q[3];
sx q[3];
rz(-1.8835246) q[3];
sx q[3];
rz(-2.6261331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];