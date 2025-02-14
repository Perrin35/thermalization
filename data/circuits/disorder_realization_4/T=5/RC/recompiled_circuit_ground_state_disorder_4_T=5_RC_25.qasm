OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0535799) q[0];
sx q[0];
rz(-0.3126643) q[0];
sx q[0];
rz(3.0076658) q[0];
rz(-0.0066095134) q[1];
sx q[1];
rz(4.7314965) q[1];
sx q[1];
rz(9.849698) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6962587) q[0];
sx q[0];
rz(-2.9888392) q[0];
sx q[0];
rz(0.85102083) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38376439) q[2];
sx q[2];
rz(-0.94049938) q[2];
sx q[2];
rz(-0.074873222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1462789) q[1];
sx q[1];
rz(-1.2321207) q[1];
sx q[1];
rz(0.92471735) q[1];
rz(-pi) q[2];
x q[2];
rz(2.22394) q[3];
sx q[3];
rz(-2.1820531) q[3];
sx q[3];
rz(-2.8070246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46140823) q[2];
sx q[2];
rz(-0.37698656) q[2];
sx q[2];
rz(3.027463) q[2];
rz(2.7528609) q[3];
sx q[3];
rz(-2.8753493) q[3];
sx q[3];
rz(1.8069327) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1330133) q[0];
sx q[0];
rz(-2.7393434) q[0];
sx q[0];
rz(-0.026799686) q[0];
rz(1.1773479) q[1];
sx q[1];
rz(-0.64759308) q[1];
sx q[1];
rz(0.037638232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5821032) q[0];
sx q[0];
rz(-2.3452098) q[0];
sx q[0];
rz(-1.5519928) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22372795) q[2];
sx q[2];
rz(-1.7729974) q[2];
sx q[2];
rz(-3.1333609) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85559713) q[1];
sx q[1];
rz(-1.0904795) q[1];
sx q[1];
rz(0.91676401) q[1];
rz(-0.82017558) q[3];
sx q[3];
rz(-1.7595152) q[3];
sx q[3];
rz(0.80126002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2521952) q[2];
sx q[2];
rz(-0.9844206) q[2];
sx q[2];
rz(-2.9241015) q[2];
rz(2.7018231) q[3];
sx q[3];
rz(-1.7871126) q[3];
sx q[3];
rz(-0.10194889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52461034) q[0];
sx q[0];
rz(-0.007402448) q[0];
sx q[0];
rz(-2.5819085) q[0];
rz(-0.89645487) q[1];
sx q[1];
rz(-2.6694916) q[1];
sx q[1];
rz(2.7339973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64227885) q[0];
sx q[0];
rz(-0.98897213) q[0];
sx q[0];
rz(0.13811843) q[0];
rz(0.90089519) q[2];
sx q[2];
rz(-1.5172361) q[2];
sx q[2];
rz(0.49240193) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36440186) q[1];
sx q[1];
rz(-2.6620416) q[1];
sx q[1];
rz(-0.24141356) q[1];
x q[2];
rz(0.30826195) q[3];
sx q[3];
rz(-2.097766) q[3];
sx q[3];
rz(1.4013724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.018230351) q[2];
sx q[2];
rz(-1.9841649) q[2];
sx q[2];
rz(1.0250214) q[2];
rz(-1.4347264) q[3];
sx q[3];
rz(-2.6985109) q[3];
sx q[3];
rz(-2.4923435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0593112) q[0];
sx q[0];
rz(-0.27114961) q[0];
sx q[0];
rz(0.49736381) q[0];
rz(1.2936032) q[1];
sx q[1];
rz(-1.0062287) q[1];
sx q[1];
rz(1.9088378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13834223) q[0];
sx q[0];
rz(-1.7742533) q[0];
sx q[0];
rz(0.5749216) q[0];
rz(-0.7035162) q[2];
sx q[2];
rz(-1.8986964) q[2];
sx q[2];
rz(0.10645535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1453443) q[1];
sx q[1];
rz(-1.323902) q[1];
sx q[1];
rz(1.3614142) q[1];
rz(-pi) q[2];
rz(-2.74128) q[3];
sx q[3];
rz(-1.5535206) q[3];
sx q[3];
rz(-0.13875419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5645912) q[2];
sx q[2];
rz(-2.5265054) q[2];
sx q[2];
rz(2.9992529) q[2];
rz(-1.9650991) q[3];
sx q[3];
rz(-1.0005955) q[3];
sx q[3];
rz(-2.71079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7560526) q[0];
sx q[0];
rz(-2.1591594) q[0];
sx q[0];
rz(2.8681712) q[0];
rz(-2.7418819) q[1];
sx q[1];
rz(-2.3276261) q[1];
sx q[1];
rz(-2.8846557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169835) q[0];
sx q[0];
rz(-2.5053609) q[0];
sx q[0];
rz(-1.0093635) q[0];
x q[1];
rz(-1.1838205) q[2];
sx q[2];
rz(-1.8024901) q[2];
sx q[2];
rz(0.53342067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94019425) q[1];
sx q[1];
rz(-1.640854) q[1];
sx q[1];
rz(0.045658535) q[1];
rz(-0.89135783) q[3];
sx q[3];
rz(-1.3472424) q[3];
sx q[3];
rz(2.4628061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5520681) q[2];
sx q[2];
rz(-2.7858211) q[2];
sx q[2];
rz(2.3815928) q[2];
rz(0.99203569) q[3];
sx q[3];
rz(-1.2090679) q[3];
sx q[3];
rz(-1.1442643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29042596) q[0];
sx q[0];
rz(-2.5100584) q[0];
sx q[0];
rz(0.60612154) q[0];
rz(-1.0468613) q[1];
sx q[1];
rz(-1.4907587) q[1];
sx q[1];
rz(-3.0643588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4299049) q[0];
sx q[0];
rz(-0.042379286) q[0];
sx q[0];
rz(1.2689356) q[0];
rz(-pi) q[1];
rz(1.0995819) q[2];
sx q[2];
rz(-1.856664) q[2];
sx q[2];
rz(0.8980823) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.162335) q[1];
sx q[1];
rz(-1.1282776) q[1];
sx q[1];
rz(2.6975432) q[1];
rz(0.66291727) q[3];
sx q[3];
rz(-1.4203216) q[3];
sx q[3];
rz(1.6126954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91654009) q[2];
sx q[2];
rz(-2.4317135) q[2];
sx q[2];
rz(-0.27099657) q[2];
rz(0.58800507) q[3];
sx q[3];
rz(-2.3376412) q[3];
sx q[3];
rz(-0.40905455) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005017) q[0];
sx q[0];
rz(-1.5301457) q[0];
sx q[0];
rz(-2.8588168) q[0];
rz(-2.458789) q[1];
sx q[1];
rz(-0.86137259) q[1];
sx q[1];
rz(2.2883794) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.565292) q[0];
sx q[0];
rz(-1.9423663) q[0];
sx q[0];
rz(2.316733) q[0];
rz(-pi) q[1];
rz(-0.81099895) q[2];
sx q[2];
rz(-1.930522) q[2];
sx q[2];
rz(-0.0051509858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4567199) q[1];
sx q[1];
rz(-1.3383075) q[1];
sx q[1];
rz(2.9732735) q[1];
x q[2];
rz(0.12577943) q[3];
sx q[3];
rz(-1.7656293) q[3];
sx q[3];
rz(1.3013713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0679438) q[2];
sx q[2];
rz(-2.6239519) q[2];
sx q[2];
rz(-2.850387) q[2];
rz(1.9368885) q[3];
sx q[3];
rz(-0.57345814) q[3];
sx q[3];
rz(1.5709491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35029995) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(-0.6268025) q[0];
rz(1.0621915) q[1];
sx q[1];
rz(-2.3419582) q[1];
sx q[1];
rz(2.3041384) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8663794) q[0];
sx q[0];
rz(-1.4555706) q[0];
sx q[0];
rz(2.000745) q[0];
rz(2.8026514) q[2];
sx q[2];
rz(-1.3474479) q[2];
sx q[2];
rz(-1.3517584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9161842) q[1];
sx q[1];
rz(-1.7042158) q[1];
sx q[1];
rz(-1.0561942) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4222048) q[3];
sx q[3];
rz(-2.0839543) q[3];
sx q[3];
rz(-3.0363014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8219139) q[2];
sx q[2];
rz(-1.9210812) q[2];
sx q[2];
rz(-1.0411881) q[2];
rz(0.90211165) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(2.5743217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64004552) q[0];
sx q[0];
rz(-2.8431659) q[0];
sx q[0];
rz(-0.10064594) q[0];
rz(1.3238662) q[1];
sx q[1];
rz(-0.83815014) q[1];
sx q[1];
rz(-0.59555882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5685265) q[0];
sx q[0];
rz(-1.0254481) q[0];
sx q[0];
rz(-0.57485707) q[0];
x q[1];
rz(-0.20636348) q[2];
sx q[2];
rz(-2.2225755) q[2];
sx q[2];
rz(-0.18557063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7157195) q[1];
sx q[1];
rz(-1.3242448) q[1];
sx q[1];
rz(2.4049525) q[1];
x q[2];
rz(2.3378701) q[3];
sx q[3];
rz(-1.6801843) q[3];
sx q[3];
rz(0.74679971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9969534) q[2];
sx q[2];
rz(-2.3914631) q[2];
sx q[2];
rz(1.7993125) q[2];
rz(-2.4053549) q[3];
sx q[3];
rz(-2.8056371) q[3];
sx q[3];
rz(3.0430702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017224273) q[0];
sx q[0];
rz(-0.65171826) q[0];
sx q[0];
rz(2.8454054) q[0];
rz(-1.724297) q[1];
sx q[1];
rz(-0.92767757) q[1];
sx q[1];
rz(2.9790624) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70746391) q[0];
sx q[0];
rz(-1.59701) q[0];
sx q[0];
rz(-3.1349584) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1549781) q[2];
sx q[2];
rz(-2.5454945) q[2];
sx q[2];
rz(1.7510027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.988058) q[1];
sx q[1];
rz(-2.9079707) q[1];
sx q[1];
rz(-2.5070666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3188884) q[3];
sx q[3];
rz(-0.71559042) q[3];
sx q[3];
rz(2.0377318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6235003) q[2];
sx q[2];
rz(-1.9237498) q[2];
sx q[2];
rz(0.45563844) q[2];
rz(-1.0738922) q[3];
sx q[3];
rz(-2.5116601) q[3];
sx q[3];
rz(0.58551252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475273) q[0];
sx q[0];
rz(-1.8210664) q[0];
sx q[0];
rz(-0.6022712) q[0];
rz(-2.8216254) q[1];
sx q[1];
rz(-1.1155557) q[1];
sx q[1];
rz(-1.3394578) q[1];
rz(-2.9721369) q[2];
sx q[2];
rz(-2.8002938) q[2];
sx q[2];
rz(2.4677966) q[2];
rz(0.77508853) q[3];
sx q[3];
rz(-1.6476484) q[3];
sx q[3];
rz(-2.021029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
