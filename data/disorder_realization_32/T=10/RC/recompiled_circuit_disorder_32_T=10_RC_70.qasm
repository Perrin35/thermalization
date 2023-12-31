OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30713233) q[0];
sx q[0];
rz(-1.6856598) q[0];
sx q[0];
rz(-2.1644724) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7726937) q[2];
sx q[2];
rz(-1.2374094) q[2];
sx q[2];
rz(2.8436529) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6160994) q[1];
sx q[1];
rz(-1.6615189) q[1];
sx q[1];
rz(-0.70848042) q[1];
x q[2];
rz(1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(-1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(2.853945) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-0.57759181) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47050414) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(-1.6862306) q[0];
x q[1];
rz(1.8446484) q[2];
sx q[2];
rz(-1.2784064) q[2];
sx q[2];
rz(-1.3016303) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2240552) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(0.097150306) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9868449) q[3];
sx q[3];
rz(-1.4156716) q[3];
sx q[3];
rz(-2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(1.0822901) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.9972237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9632918) q[0];
sx q[0];
rz(-0.68094567) q[0];
sx q[0];
rz(2.2292024) q[0];
rz(-pi) q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(-0.70659107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5971165) q[1];
sx q[1];
rz(-2.5529478) q[1];
sx q[1];
rz(-2.4878923) q[1];
rz(-pi) q[2];
rz(-1.5677489) q[3];
sx q[3];
rz(-1.8584195) q[3];
sx q[3];
rz(-2.0730719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3545761) q[0];
sx q[0];
rz(-2.0873318) q[0];
sx q[0];
rz(1.2660962) q[0];
rz(-0.86187141) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(2.9574403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0419473) q[1];
sx q[1];
rz(-1.6640267) q[1];
sx q[1];
rz(-0.16361841) q[1];
rz(-pi) q[2];
rz(-0.99289258) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(-2.2955017) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-2.3748421) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2131166) q[0];
sx q[0];
rz(-2.1931744) q[0];
sx q[0];
rz(1.9268131) q[0];
x q[1];
rz(-0.93047662) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(-0.1375246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5366718) q[1];
sx q[1];
rz(-0.39502883) q[1];
sx q[1];
rz(-3.1035963) q[1];
rz(-pi) q[2];
rz(-2.1122123) q[3];
sx q[3];
rz(-2.5684528) q[3];
sx q[3];
rz(-1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(1.0995964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.106364) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(2.7307672) q[0];
x q[1];
rz(-1.2560558) q[2];
sx q[2];
rz(-2.5883) q[2];
sx q[2];
rz(1.7994583) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9601599) q[1];
sx q[1];
rz(-0.90248855) q[1];
sx q[1];
rz(1.798435) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.908329) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(1.69676) q[0];
rz(-pi) q[1];
rz(-1.332875) q[2];
sx q[2];
rz(-0.750713) q[2];
sx q[2];
rz(-2.6964292) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8017756) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(2.9031309) q[1];
rz(-1.8623452) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(-0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(2.3106993) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6458994) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(-0.010852531) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56284833) q[0];
sx q[0];
rz(-1.0283854) q[0];
sx q[0];
rz(1.5777274) q[0];
rz(-pi) q[1];
x q[1];
rz(1.79004) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(3.0041681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9353232) q[1];
sx q[1];
rz(-0.85201293) q[1];
sx q[1];
rz(-2.1410336) q[1];
rz(-pi) q[2];
rz(1.9803489) q[3];
sx q[3];
rz(-1.1405986) q[3];
sx q[3];
rz(-0.48914117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(-0.44257277) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.7907422) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.6500236) q[0];
sx q[0];
rz(-1.4863192) q[0];
rz(-1.4883792) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(-1.272162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1119712) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(3.1236468) q[1];
rz(-pi) q[2];
rz(-1.7926932) q[3];
sx q[3];
rz(-1.4699234) q[3];
sx q[3];
rz(-1.5125121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1039625) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(-0.9141586) q[0];
x q[1];
rz(3.0868297) q[2];
sx q[2];
rz(-2.7191396) q[2];
sx q[2];
rz(-1.7516608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1093724) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(-2.1857775) q[1];
x q[2];
rz(1.1630837) q[3];
sx q[3];
rz(-2.793503) q[3];
sx q[3];
rz(-0.40504328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(-2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8511843) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-1.2693263) q[2];
sx q[2];
rz(-1.3133247) q[2];
sx q[2];
rz(1.7016344) q[2];
rz(-0.11733304) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
