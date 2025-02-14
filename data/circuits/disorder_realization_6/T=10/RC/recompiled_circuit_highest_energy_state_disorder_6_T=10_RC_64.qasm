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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(0.52752703) q[1];
sx q[1];
rz(-2.1105284) q[1];
sx q[1];
rz(0.48479015) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82085375) q[0];
sx q[0];
rz(-1.0620097) q[0];
sx q[0];
rz(0.15055045) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3622857) q[2];
sx q[2];
rz(-1.2504842) q[2];
sx q[2];
rz(-0.28398289) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97580662) q[1];
sx q[1];
rz(-0.96250223) q[1];
sx q[1];
rz(-0.071578524) q[1];
rz(-1.6438585) q[3];
sx q[3];
rz(-1.2803852) q[3];
sx q[3];
rz(-0.452347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63554865) q[2];
sx q[2];
rz(-0.94749331) q[2];
sx q[2];
rz(-0.49723899) q[2];
rz(-0.55073589) q[3];
sx q[3];
rz(-0.68998718) q[3];
sx q[3];
rz(-1.1373038) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15744844) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(-0.76341158) q[0];
rz(-2.5570671) q[1];
sx q[1];
rz(-1.3817363) q[1];
sx q[1];
rz(1.0113299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4268734) q[0];
sx q[0];
rz(-1.4397286) q[0];
sx q[0];
rz(3.0253973) q[0];
rz(-pi) q[1];
rz(-2.6945128) q[2];
sx q[2];
rz(-2.8805388) q[2];
sx q[2];
rz(-2.3787468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3015131) q[1];
sx q[1];
rz(-0.40427819) q[1];
sx q[1];
rz(2.1188645) q[1];
rz(0.92722757) q[3];
sx q[3];
rz(-2.5303741) q[3];
sx q[3];
rz(-0.54655308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4804907) q[2];
sx q[2];
rz(-1.9384408) q[2];
sx q[2];
rz(-1.0002381) q[2];
rz(-0.58146042) q[3];
sx q[3];
rz(-0.73271078) q[3];
sx q[3];
rz(1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8345399) q[0];
sx q[0];
rz(-0.59891278) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(0.3071951) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(0.44816005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7355928) q[0];
sx q[0];
rz(-1.4383398) q[0];
sx q[0];
rz(0.51901061) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9314172) q[2];
sx q[2];
rz(-0.54979392) q[2];
sx q[2];
rz(0.87042337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8413755) q[1];
sx q[1];
rz(-0.63134495) q[1];
sx q[1];
rz(-1.6168827) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5968733) q[3];
sx q[3];
rz(-1.136565) q[3];
sx q[3];
rz(-2.6573879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8223411) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(0.59784293) q[2];
rz(2.6939825) q[3];
sx q[3];
rz(-1.8748583) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.7287801) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(-0.70984167) q[0];
rz(2.1624883) q[1];
sx q[1];
rz(-2.282228) q[1];
sx q[1];
rz(1.8202579) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431491) q[0];
sx q[0];
rz(-1.51529) q[0];
sx q[0];
rz(2.4401791) q[0];
rz(-pi) q[1];
rz(0.780402) q[2];
sx q[2];
rz(-1.4596887) q[2];
sx q[2];
rz(-1.358913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37934694) q[1];
sx q[1];
rz(-1.3031465) q[1];
sx q[1];
rz(0.91581151) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31339733) q[3];
sx q[3];
rz(-1.0492696) q[3];
sx q[3];
rz(2.1300621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9909624) q[2];
sx q[2];
rz(-1.7143098) q[2];
sx q[2];
rz(2.4347351) q[2];
rz(-1.3458716) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(2.4367387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83955806) q[0];
sx q[0];
rz(-2.3607881) q[0];
sx q[0];
rz(1.5020405) q[0];
rz(1.543965) q[1];
sx q[1];
rz(-0.85999703) q[1];
sx q[1];
rz(1.647321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92163819) q[0];
sx q[0];
rz(-2.1412379) q[0];
sx q[0];
rz(-0.0019689671) q[0];
rz(1.5780894) q[2];
sx q[2];
rz(-0.97541684) q[2];
sx q[2];
rz(-3.1356168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8710821) q[1];
sx q[1];
rz(-1.8746261) q[1];
sx q[1];
rz(-2.8302776) q[1];
rz(0.95559883) q[3];
sx q[3];
rz(-1.5892059) q[3];
sx q[3];
rz(1.9083202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8813701) q[2];
sx q[2];
rz(-1.8689195) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(-2.6327366) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(-1.2671027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.59638554) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(2.9811356) q[0];
rz(2.4682553) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(2.014726) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.883553) q[0];
sx q[0];
rz(-0.48012381) q[0];
sx q[0];
rz(-2.2715946) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2484064) q[2];
sx q[2];
rz(-2.3087072) q[2];
sx q[2];
rz(-2.5423342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8452478) q[1];
sx q[1];
rz(-2.5242477) q[1];
sx q[1];
rz(-1.9250573) q[1];
x q[2];
rz(1.2282727) q[3];
sx q[3];
rz(-0.5687826) q[3];
sx q[3];
rz(2.3186213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0141853) q[2];
sx q[2];
rz(-1.5994453) q[2];
sx q[2];
rz(2.3678153) q[2];
rz(2.1524147) q[3];
sx q[3];
rz(-2.7343605) q[3];
sx q[3];
rz(2.0686843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147375) q[0];
sx q[0];
rz(-1.4829153) q[0];
sx q[0];
rz(-0.7731272) q[0];
rz(1.5559366) q[1];
sx q[1];
rz(-1.2890041) q[1];
sx q[1];
rz(-1.1770491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0814514) q[0];
sx q[0];
rz(-0.8340652) q[0];
sx q[0];
rz(0.042515083) q[0];
rz(-pi) q[1];
rz(-0.99780166) q[2];
sx q[2];
rz(-1.3062451) q[2];
sx q[2];
rz(-0.71392347) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.24304535) q[1];
sx q[1];
rz(-2.4345401) q[1];
sx q[1];
rz(-1.1416332) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7103393) q[3];
sx q[3];
rz(-1.7467032) q[3];
sx q[3];
rz(-1.6493451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3561463) q[2];
sx q[2];
rz(-0.36474228) q[2];
sx q[2];
rz(2.6981603) q[2];
rz(-0.38608471) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(-0.20080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3767553) q[0];
sx q[0];
rz(-1.5926462) q[0];
sx q[0];
rz(0.03373294) q[0];
rz(0.69817606) q[1];
sx q[1];
rz(-2.5444784) q[1];
sx q[1];
rz(-2.4765292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2584943) q[0];
sx q[0];
rz(-1.8332714) q[0];
sx q[0];
rz(-1.1061944) q[0];
x q[1];
rz(2.0940789) q[2];
sx q[2];
rz(-2.3436597) q[2];
sx q[2];
rz(2.8454859) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60282502) q[1];
sx q[1];
rz(-2.1884349) q[1];
sx q[1];
rz(-2.5185939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5916717) q[3];
sx q[3];
rz(-1.2513147) q[3];
sx q[3];
rz(3.0221981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.72208059) q[2];
sx q[2];
rz(-2.2245912) q[2];
sx q[2];
rz(-1.7783995) q[2];
rz(2.6204387) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9687965) q[0];
sx q[0];
rz(-1.7954614) q[0];
sx q[0];
rz(-2.3222493) q[0];
rz(-0.68483812) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(3.0192764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8259895) q[0];
sx q[0];
rz(-1.9902744) q[0];
sx q[0];
rz(-0.88780888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38057957) q[2];
sx q[2];
rz(-2.1316445) q[2];
sx q[2];
rz(-0.83319908) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.015339) q[1];
sx q[1];
rz(-0.97927927) q[1];
sx q[1];
rz(-1.1777608) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4966481) q[3];
sx q[3];
rz(-2.44584) q[3];
sx q[3];
rz(-0.49518817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(0.86755794) q[2];
rz(-1.3982754) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(2.5625663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45359465) q[0];
sx q[0];
rz(-1.9168251) q[0];
sx q[0];
rz(2.0523742) q[0];
rz(0.047529686) q[1];
sx q[1];
rz(-1.0462953) q[1];
sx q[1];
rz(0.68738031) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9935038) q[0];
sx q[0];
rz(-1.5951135) q[0];
sx q[0];
rz(2.0612102) q[0];
x q[1];
rz(2.8609852) q[2];
sx q[2];
rz(-1.1550552) q[2];
sx q[2];
rz(-3.0688697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78391548) q[1];
sx q[1];
rz(-1.1184753) q[1];
sx q[1];
rz(0.6953101) q[1];
rz(-pi) q[2];
rz(-0.093133868) q[3];
sx q[3];
rz(-2.2054005) q[3];
sx q[3];
rz(-0.30090162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99531949) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(1.6590365) q[2];
rz(2.4847374) q[3];
sx q[3];
rz(-1.992179) q[3];
sx q[3];
rz(0.38203865) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.252608) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(2.5629015) q[1];
sx q[1];
rz(-0.83732579) q[1];
sx q[1];
rz(2.8566828) q[1];
rz(-1.4877697) q[2];
sx q[2];
rz(-1.2157233) q[2];
sx q[2];
rz(0.12479648) q[2];
rz(-1.4887078) q[3];
sx q[3];
rz(-1.0552867) q[3];
sx q[3];
rz(-1.1724006) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
