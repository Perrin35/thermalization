OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7772726) q[0];
sx q[0];
rz(-0.84982818) q[0];
sx q[0];
rz(-1.538869) q[0];
rz(2.377805) q[1];
sx q[1];
rz(-3.0883767) q[1];
sx q[1];
rz(2.3071212) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6042701) q[0];
sx q[0];
rz(-1.5198845) q[0];
sx q[0];
rz(2.0432274) q[0];
rz(2.3896289) q[2];
sx q[2];
rz(-2.3533216) q[2];
sx q[2];
rz(1.0984818) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.758513) q[1];
sx q[1];
rz(-2.5316409) q[1];
sx q[1];
rz(-1.3332696) q[1];
rz(-pi) q[2];
rz(-2.2356307) q[3];
sx q[3];
rz(-0.49091456) q[3];
sx q[3];
rz(-2.5708446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6795464) q[2];
sx q[2];
rz(-1.5507689) q[2];
sx q[2];
rz(0.66590434) q[2];
rz(0.42689651) q[3];
sx q[3];
rz(-2.175943) q[3];
sx q[3];
rz(0.20195937) q[3];
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
rz(2.9694328) q[0];
sx q[0];
rz(-2.1803624) q[0];
sx q[0];
rz(0.52312553) q[0];
rz(0.41921774) q[1];
sx q[1];
rz(-2.9328177) q[1];
sx q[1];
rz(-2.8255393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284477) q[0];
sx q[0];
rz(-2.4372167) q[0];
sx q[0];
rz(2.8718642) q[0];
x q[1];
rz(0.14618539) q[2];
sx q[2];
rz(-1.5461032) q[2];
sx q[2];
rz(3.0865412) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.7853941) q[1];
sx q[1];
rz(-1.0266582) q[1];
sx q[1];
rz(-0.35299976) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7394946) q[3];
sx q[3];
rz(-1.1505812) q[3];
sx q[3];
rz(1.4762154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6009964) q[2];
sx q[2];
rz(-0.54800802) q[2];
sx q[2];
rz(-2.0151286) q[2];
rz(-0.90538853) q[3];
sx q[3];
rz(-1.0631881) q[3];
sx q[3];
rz(-1.5409013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044540502) q[0];
sx q[0];
rz(-1.4690761) q[0];
sx q[0];
rz(1.4918406) q[0];
rz(0.39332321) q[1];
sx q[1];
rz(-1.2231188) q[1];
sx q[1];
rz(-0.16600674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1309083) q[0];
sx q[0];
rz(-1.8839399) q[0];
sx q[0];
rz(-0.74096276) q[0];
rz(1.7766062) q[2];
sx q[2];
rz(-2.4327288) q[2];
sx q[2];
rz(1.0212932) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5854726) q[1];
sx q[1];
rz(-2.2801546) q[1];
sx q[1];
rz(1.8147537) q[1];
x q[2];
rz(2.7083206) q[3];
sx q[3];
rz(-1.6973933) q[3];
sx q[3];
rz(-1.8357852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7321221) q[2];
sx q[2];
rz(-2.8889443) q[2];
sx q[2];
rz(2.9518993) q[2];
rz(0.41483748) q[3];
sx q[3];
rz(-1.8388351) q[3];
sx q[3];
rz(-0.10492575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6787427) q[0];
sx q[0];
rz(-0.78063709) q[0];
sx q[0];
rz(-1.8829874) q[0];
rz(1.7557834) q[1];
sx q[1];
rz(-0.36711991) q[1];
sx q[1];
rz(0.78305903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5393128) q[0];
sx q[0];
rz(-3.1010744) q[0];
sx q[0];
rz(0.65538533) q[0];
rz(-1.5668529) q[2];
sx q[2];
rz(-0.76216216) q[2];
sx q[2];
rz(0.54715014) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9077346) q[1];
sx q[1];
rz(-2.0480262) q[1];
sx q[1];
rz(-2.1388172) q[1];
rz(-1.542605) q[3];
sx q[3];
rz(-1.2148093) q[3];
sx q[3];
rz(0.062075768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2087525) q[2];
sx q[2];
rz(-1.1609266) q[2];
sx q[2];
rz(2.0483542) q[2];
rz(2.7161963) q[3];
sx q[3];
rz(-1.930611) q[3];
sx q[3];
rz(-1.0824599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071674034) q[0];
sx q[0];
rz(-2.2921083) q[0];
sx q[0];
rz(1.7639532) q[0];
rz(-2.2397974) q[1];
sx q[1];
rz(-1.8865562) q[1];
sx q[1];
rz(-0.76595438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655584) q[0];
sx q[0];
rz(-1.9355825) q[0];
sx q[0];
rz(0.64393142) q[0];
rz(-2.7904205) q[2];
sx q[2];
rz(-2.4173556) q[2];
sx q[2];
rz(-1.8623095) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53772012) q[1];
sx q[1];
rz(-2.0687547) q[1];
sx q[1];
rz(1.3977951) q[1];
rz(1.1328825) q[3];
sx q[3];
rz(-1.8275598) q[3];
sx q[3];
rz(-1.6516466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9506266) q[2];
sx q[2];
rz(-0.85056716) q[2];
sx q[2];
rz(-1.4858049) q[2];
rz(-2.0048257) q[3];
sx q[3];
rz(-0.4509238) q[3];
sx q[3];
rz(2.3740812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9039827) q[0];
sx q[0];
rz(-0.5466277) q[0];
sx q[0];
rz(2.9055063) q[0];
rz(-1.4769332) q[1];
sx q[1];
rz(-1.4918985) q[1];
sx q[1];
rz(-1.9794855) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6593137) q[0];
sx q[0];
rz(-1.3516607) q[0];
sx q[0];
rz(-2.3740039) q[0];
rz(-pi) q[1];
rz(0.40923758) q[2];
sx q[2];
rz(-1.6820568) q[2];
sx q[2];
rz(-1.2986599) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2233189) q[1];
sx q[1];
rz(-1.9293203) q[1];
sx q[1];
rz(-2.8467092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.052836494) q[3];
sx q[3];
rz(-1.5457258) q[3];
sx q[3];
rz(-2.3468338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7621925) q[2];
sx q[2];
rz(-0.47319943) q[2];
sx q[2];
rz(-0.0089448373) q[2];
rz(1.829932) q[3];
sx q[3];
rz(-1.0511755) q[3];
sx q[3];
rz(-1.293762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0084429) q[0];
sx q[0];
rz(-2.382998) q[0];
sx q[0];
rz(-1.5822509) q[0];
rz(2.3186231) q[1];
sx q[1];
rz(-0.4987078) q[1];
sx q[1];
rz(-1.4272326) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5873272) q[0];
sx q[0];
rz(-1.1236941) q[0];
sx q[0];
rz(-2.7914911) q[0];
rz(-0.64325579) q[2];
sx q[2];
rz(-1.3998886) q[2];
sx q[2];
rz(2.3115445) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8135951) q[1];
sx q[1];
rz(-2.5668721) q[1];
sx q[1];
rz(-1.5411472) q[1];
x q[2];
rz(1.6839833) q[3];
sx q[3];
rz(-1.9838399) q[3];
sx q[3];
rz(2.7922009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83671776) q[2];
sx q[2];
rz(-1.0739526) q[2];
sx q[2];
rz(2.7433336) q[2];
rz(-2.7598329) q[3];
sx q[3];
rz(-1.7002117) q[3];
sx q[3];
rz(2.1329342) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1036103) q[0];
sx q[0];
rz(-1.4382265) q[0];
sx q[0];
rz(1.2305228) q[0];
rz(-1.1331753) q[1];
sx q[1];
rz(-0.76825348) q[1];
sx q[1];
rz(-0.59655985) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2418848) q[0];
sx q[0];
rz(-1.5800137) q[0];
sx q[0];
rz(-0.25141821) q[0];
rz(-pi) q[1];
rz(1.2136772) q[2];
sx q[2];
rz(-2.6304863) q[2];
sx q[2];
rz(2.7765283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.34932) q[1];
sx q[1];
rz(-1.8304003) q[1];
sx q[1];
rz(-2.1672897) q[1];
rz(0.61847134) q[3];
sx q[3];
rz(-1.2534646) q[3];
sx q[3];
rz(-2.0813074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29264984) q[2];
sx q[2];
rz(-0.611648) q[2];
sx q[2];
rz(1.6742279) q[2];
rz(-2.8749706) q[3];
sx q[3];
rz(-2.1610503) q[3];
sx q[3];
rz(-2.9149741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6999231) q[0];
sx q[0];
rz(-3.141098) q[0];
sx q[0];
rz(2.0200404) q[0];
rz(-0.17555217) q[1];
sx q[1];
rz(-2.4708864) q[1];
sx q[1];
rz(0.49351722) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8275857) q[0];
sx q[0];
rz(-0.092484154) q[0];
sx q[0];
rz(0.88351639) q[0];
rz(-2.001112) q[2];
sx q[2];
rz(-2.3386933) q[2];
sx q[2];
rz(-1.0969298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35084823) q[1];
sx q[1];
rz(-2.1022294) q[1];
sx q[1];
rz(2.8709251) q[1];
rz(-pi) q[2];
rz(1.4021691) q[3];
sx q[3];
rz(-2.6624749) q[3];
sx q[3];
rz(-1.7134596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4404122) q[2];
sx q[2];
rz(-2.7737325) q[2];
sx q[2];
rz(2.5555447) q[2];
rz(1.690257) q[3];
sx q[3];
rz(-1.3380932) q[3];
sx q[3];
rz(0.51226789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.119809) q[0];
sx q[0];
rz(-1.021294) q[0];
sx q[0];
rz(0.16948608) q[0];
rz(0.64300621) q[1];
sx q[1];
rz(-0.87130672) q[1];
sx q[1];
rz(-0.51173425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62637037) q[0];
sx q[0];
rz(-3.0128509) q[0];
sx q[0];
rz(-2.4730352) q[0];
rz(-pi) q[1];
rz(-1.636347) q[2];
sx q[2];
rz(-1.3093072) q[2];
sx q[2];
rz(-2.477867) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4525254) q[1];
sx q[1];
rz(-2.2241909) q[1];
sx q[1];
rz(0.19028581) q[1];
rz(0.2318895) q[3];
sx q[3];
rz(-2.0733193) q[3];
sx q[3];
rz(-2.1351867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6181651) q[2];
sx q[2];
rz(-0.7343505) q[2];
sx q[2];
rz(0.0283453) q[2];
rz(-0.47447765) q[3];
sx q[3];
rz(-2.7826383) q[3];
sx q[3];
rz(-1.0396144) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9356553) q[0];
sx q[0];
rz(-1.7556998) q[0];
sx q[0];
rz(2.6082267) q[0];
rz(-0.43279303) q[1];
sx q[1];
rz(-1.5569893) q[1];
sx q[1];
rz(2.5098262) q[1];
rz(1.648394) q[2];
sx q[2];
rz(-1.785555) q[2];
sx q[2];
rz(2.6746489) q[2];
rz(-1.2148874) q[3];
sx q[3];
rz(-2.9199745) q[3];
sx q[3];
rz(-1.5501116) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
