OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18093827) q[0];
sx q[0];
rz(-3.0783983) q[0];
sx q[0];
rz(0.59046459) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(1.7658748) q[1];
sx q[1];
rz(9.4212846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1549415) q[0];
sx q[0];
rz(-1.1677603) q[0];
sx q[0];
rz(2.0986845) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3071877) q[2];
sx q[2];
rz(-2.2279539) q[2];
sx q[2];
rz(3.1176569) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3319015) q[1];
sx q[1];
rz(-1.594048) q[1];
sx q[1];
rz(-1.5493259) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6925647) q[3];
sx q[3];
rz(-2.1204344) q[3];
sx q[3];
rz(-2.6688936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8840088) q[2];
sx q[2];
rz(-2.282892) q[2];
sx q[2];
rz(1.0165455) q[2];
rz(0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(-1.6421002) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4676374) q[0];
sx q[0];
rz(-1.6683945) q[0];
sx q[0];
rz(3.0067645) q[0];
rz(0.22678953) q[1];
sx q[1];
rz(-3.0786381) q[1];
sx q[1];
rz(2.8917868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9591208) q[0];
sx q[0];
rz(-1.153649) q[0];
sx q[0];
rz(2.4802101) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.05257865) q[2];
sx q[2];
rz(-2.0705531) q[2];
sx q[2];
rz(1.5165968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0909101) q[1];
sx q[1];
rz(-1.5600648) q[1];
sx q[1];
rz(-3.0640567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2050912) q[3];
sx q[3];
rz(-1.6964751) q[3];
sx q[3];
rz(-0.92363908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9164385) q[2];
sx q[2];
rz(-0.068000451) q[2];
sx q[2];
rz(2.1615243) q[2];
rz(0.89503908) q[3];
sx q[3];
rz(-2.3990302) q[3];
sx q[3];
rz(1.9612954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5535683) q[0];
sx q[0];
rz(-2.6346485) q[0];
sx q[0];
rz(-1.5356327) q[0];
rz(-1.5060679) q[1];
sx q[1];
rz(-2.3160544) q[1];
sx q[1];
rz(-0.95605409) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0515389) q[0];
sx q[0];
rz(-0.23583007) q[0];
sx q[0];
rz(-0.43311849) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73101513) q[2];
sx q[2];
rz(-1.6201978) q[2];
sx q[2];
rz(1.4412944) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29143476) q[1];
sx q[1];
rz(-2.4031284) q[1];
sx q[1];
rz(1.4126331) q[1];
rz(-pi) q[2];
rz(0.51436971) q[3];
sx q[3];
rz(-1.8764495) q[3];
sx q[3];
rz(-2.5955868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.306159) q[2];
sx q[2];
rz(-2.4296438) q[2];
sx q[2];
rz(-0.63208675) q[2];
rz(2.2575374) q[3];
sx q[3];
rz(-3.1243117) q[3];
sx q[3];
rz(0.89068762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4759091) q[0];
sx q[0];
rz(-1.8432239) q[0];
sx q[0];
rz(-0.16715288) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-2.1985168) q[1];
sx q[1];
rz(1.4080338) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65850019) q[0];
sx q[0];
rz(-1.2804872) q[0];
sx q[0];
rz(0.4418375) q[0];
rz(-pi) q[1];
rz(1.3590668) q[2];
sx q[2];
rz(-1.5671726) q[2];
sx q[2];
rz(-0.45932367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38630292) q[1];
sx q[1];
rz(-1.3805256) q[1];
sx q[1];
rz(-2.931837) q[1];
x q[2];
rz(-2.3711331) q[3];
sx q[3];
rz(-1.7826728) q[3];
sx q[3];
rz(0.4650863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7888646) q[2];
sx q[2];
rz(-3.126324) q[2];
sx q[2];
rz(-0.35066476) q[2];
rz(-0.26385012) q[3];
sx q[3];
rz(-3.1407472) q[3];
sx q[3];
rz(1.9381757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26461399) q[0];
sx q[0];
rz(-2.3961841) q[0];
sx q[0];
rz(1.718148) q[0];
rz(-2.8599332) q[1];
sx q[1];
rz(-1.5080844) q[1];
sx q[1];
rz(-1.3196094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3521101) q[0];
sx q[0];
rz(-0.69449678) q[0];
sx q[0];
rz(-0.8775893) q[0];
x q[1];
rz(0.93792589) q[2];
sx q[2];
rz(-3.1073776) q[2];
sx q[2];
rz(0.93955112) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4180243) q[1];
sx q[1];
rz(-0.71360525) q[1];
sx q[1];
rz(2.6221971) q[1];
rz(-pi) q[2];
rz(0.67992391) q[3];
sx q[3];
rz(-2.6377684) q[3];
sx q[3];
rz(-0.50407223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8772956) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(1.7128672) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-3.0890586) q[3];
sx q[3];
rz(-3.0534237) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0577211) q[0];
sx q[0];
rz(-0.11744048) q[0];
sx q[0];
rz(2.591326) q[0];
rz(2.8394207) q[1];
sx q[1];
rz(-1.7487532) q[1];
sx q[1];
rz(0.43513939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0939729) q[0];
sx q[0];
rz(-0.21417266) q[0];
sx q[0];
rz(-0.16475494) q[0];
rz(-pi) q[1];
rz(1.5717634) q[2];
sx q[2];
rz(-1.5782326) q[2];
sx q[2];
rz(-0.32583562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3966538) q[1];
sx q[1];
rz(-0.42801539) q[1];
sx q[1];
rz(-1.4465843) q[1];
rz(2.4184259) q[3];
sx q[3];
rz(-2.1373014) q[3];
sx q[3];
rz(2.1700933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2415194) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(-0.97236902) q[2];
rz(0.9995681) q[3];
sx q[3];
rz(-0.0071503706) q[3];
sx q[3];
rz(2.0991367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59357506) q[0];
sx q[0];
rz(-1.8997718) q[0];
sx q[0];
rz(-2.0752564) q[0];
rz(1.713133) q[1];
sx q[1];
rz(-2.8722873) q[1];
sx q[1];
rz(1.2880633) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748874) q[0];
sx q[0];
rz(-2.2572491) q[0];
sx q[0];
rz(1.1576325) q[0];
rz(-pi) q[1];
rz(2.3003103) q[2];
sx q[2];
rz(-1.8389353) q[2];
sx q[2];
rz(2.8644305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8180698) q[1];
sx q[1];
rz(-1.5574291) q[1];
sx q[1];
rz(-1.531015) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24434511) q[3];
sx q[3];
rz(-2.0879896) q[3];
sx q[3];
rz(-1.2822267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0224033) q[2];
sx q[2];
rz(-3.0722805) q[2];
sx q[2];
rz(-0.27528396) q[2];
rz(0.22661181) q[3];
sx q[3];
rz(-2.785502) q[3];
sx q[3];
rz(2.7039458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76448035) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(-2.5544033) q[0];
rz(-1.5234692) q[1];
sx q[1];
rz(-2.0240929) q[1];
sx q[1];
rz(-1.5706185) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29059288) q[0];
sx q[0];
rz(-2.3208566) q[0];
sx q[0];
rz(-0.59159578) q[0];
rz(-pi) q[1];
rz(-3.0745071) q[2];
sx q[2];
rz(-2.0338981) q[2];
sx q[2];
rz(2.9912134) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60698673) q[1];
sx q[1];
rz(-0.053664975) q[1];
sx q[1];
rz(2.1855658) q[1];
rz(-pi) q[2];
rz(2.0471025) q[3];
sx q[3];
rz(-0.62783754) q[3];
sx q[3];
rz(1.353457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1945343) q[2];
sx q[2];
rz(-0.18320601) q[2];
sx q[2];
rz(-1.3018695) q[2];
rz(2.7070847) q[3];
sx q[3];
rz(-0.040381581) q[3];
sx q[3];
rz(2.6935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-1.3473564) q[0];
sx q[0];
rz(-3.0257822) q[0];
sx q[0];
rz(-1.7606803) q[0];
rz(1.5877089) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(0.10996058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0123574) q[0];
sx q[0];
rz(-1.7596272) q[0];
sx q[0];
rz(-2.3055632) q[0];
rz(-pi) q[1];
rz(-2.8822363) q[2];
sx q[2];
rz(-2.4217941) q[2];
sx q[2];
rz(2.6754745) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1404071) q[1];
sx q[1];
rz(-0.81403908) q[1];
sx q[1];
rz(-3.1073529) q[1];
x q[2];
rz(-1.7305829) q[3];
sx q[3];
rz(-1.2629804) q[3];
sx q[3];
rz(0.45444835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79003698) q[2];
sx q[2];
rz(-1.0645126) q[2];
sx q[2];
rz(-0.35247394) q[2];
rz(0.6414837) q[3];
sx q[3];
rz(-3.0951169) q[3];
sx q[3];
rz(-0.36267734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8620257) q[0];
sx q[0];
rz(-2.8643705) q[0];
sx q[0];
rz(-0.56420457) q[0];
rz(-1.5203681) q[1];
sx q[1];
rz(-1.0313326) q[1];
sx q[1];
rz(3.0618111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7120651) q[0];
sx q[0];
rz(-2.2343415) q[0];
sx q[0];
rz(0.19312858) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0941102) q[2];
sx q[2];
rz(-1.6185307) q[2];
sx q[2];
rz(-1.6636696) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3120739) q[1];
sx q[1];
rz(-0.79527277) q[1];
sx q[1];
rz(-1.9327232) q[1];
rz(-pi) q[2];
rz(-2.8989927) q[3];
sx q[3];
rz(-1.4897457) q[3];
sx q[3];
rz(1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0148049) q[2];
sx q[2];
rz(-3.1316275) q[2];
sx q[2];
rz(1.0753746) q[2];
rz(0.77438313) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(0.27124852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423691) q[0];
sx q[0];
rz(-1.5726226) q[0];
sx q[0];
rz(1.5725305) q[0];
rz(2.6161999) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(1.0030307) q[2];
sx q[2];
rz(-1.7519288) q[2];
sx q[2];
rz(-1.1733423) q[2];
rz(2.5767951) q[3];
sx q[3];
rz(-2.372143) q[3];
sx q[3];
rz(-1.8134223) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
