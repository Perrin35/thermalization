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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(2.0896572) q[0];
rz(-2.9479041) q[1];
sx q[1];
rz(-2.5084578) q[1];
sx q[1];
rz(2.5685891) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6010701) q[0];
sx q[0];
rz(-2.1946476) q[0];
sx q[0];
rz(-3.0649351) q[0];
rz(-0.51197211) q[2];
sx q[2];
rz(-2.3880771) q[2];
sx q[2];
rz(1.3775502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6261648) q[1];
sx q[1];
rz(-1.0744796) q[1];
sx q[1];
rz(1.650701) q[1];
rz(0.696508) q[3];
sx q[3];
rz(-1.4918431) q[3];
sx q[3];
rz(-2.1895308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8594592) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(1.0563043) q[2];
rz(-3.1193962) q[3];
sx q[3];
rz(-2.3789417) q[3];
sx q[3];
rz(1.7215884) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2317155) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(0.12769708) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(-1.4375623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50629726) q[0];
sx q[0];
rz(-0.40182913) q[0];
sx q[0];
rz(0.70962972) q[0];
x q[1];
rz(1.5135852) q[2];
sx q[2];
rz(-1.621641) q[2];
sx q[2];
rz(-2.0432931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54894069) q[1];
sx q[1];
rz(-0.90528622) q[1];
sx q[1];
rz(1.9926461) q[1];
rz(-0.53967472) q[3];
sx q[3];
rz(-1.3840578) q[3];
sx q[3];
rz(0.58247551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-0.44542584) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-1.0990812) q[3];
sx q[3];
rz(0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55260783) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(-2.4267922) q[0];
rz(1.0572761) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(-0.32726273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875558) q[0];
sx q[0];
rz(-1.7374514) q[0];
sx q[0];
rz(2.119407) q[0];
rz(-pi) q[1];
rz(-1.0017603) q[2];
sx q[2];
rz(-1.1998917) q[2];
sx q[2];
rz(-2.6628464) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5048557) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(1.4192753) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25423519) q[3];
sx q[3];
rz(-1.2336858) q[3];
sx q[3];
rz(1.1439307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1383692) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(-1.6376015) q[2];
rz(1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(-2.159582) q[3];
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
rz(-pi) q[3];
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
rz(-1.0729436) q[0];
sx q[0];
rz(-1.2462085) q[0];
sx q[0];
rz(-0.15596998) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-0.60395423) q[1];
sx q[1];
rz(1.2329996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.975978) q[0];
sx q[0];
rz(-2.8535643) q[0];
sx q[0];
rz(2.0095129) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92278752) q[2];
sx q[2];
rz(-1.4013259) q[2];
sx q[2];
rz(3.0344935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2501331) q[1];
sx q[1];
rz(-2.1910408) q[1];
sx q[1];
rz(-0.34858443) q[1];
x q[2];
rz(-1.5560914) q[3];
sx q[3];
rz(-1.8222408) q[3];
sx q[3];
rz(2.3594643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6984581) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(-1.5114463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(-2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-0.41780892) q[0];
sx q[0];
rz(1.9885709) q[0];
rz(2.5979089) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(-2.8797454) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3735298) q[0];
sx q[0];
rz(-0.090001194) q[0];
sx q[0];
rz(-1.1885719) q[0];
rz(-pi) q[1];
rz(1.0864429) q[2];
sx q[2];
rz(-2.5005955) q[2];
sx q[2];
rz(-0.5400368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3186372) q[1];
sx q[1];
rz(-1.3483817) q[1];
sx q[1];
rz(-2.5435124) q[1];
rz(2.35942) q[3];
sx q[3];
rz(-1.7165136) q[3];
sx q[3];
rz(-0.71469939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5220945) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(-1.4052793) q[2];
rz(-0.74553472) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(-2.1982927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.140542) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(2.70128) q[1];
sx q[1];
rz(-1.3654717) q[1];
sx q[1];
rz(2.9702759) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650314) q[0];
sx q[0];
rz(-2.2269571) q[0];
sx q[0];
rz(1.879346) q[0];
x q[1];
rz(-1.2662884) q[2];
sx q[2];
rz(-2.300548) q[2];
sx q[2];
rz(-0.27308057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.379516) q[1];
sx q[1];
rz(-0.87240309) q[1];
sx q[1];
rz(-1.5980491) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8796575) q[3];
sx q[3];
rz(-0.98239693) q[3];
sx q[3];
rz(-0.0042875687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0773086) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(0.94179955) q[2];
rz(1.9866379) q[3];
sx q[3];
rz(-2.933511) q[3];
sx q[3];
rz(-1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44523859) q[0];
sx q[0];
rz(-2.0360763) q[0];
sx q[0];
rz(-0.0048986991) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-2.4536224) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94243455) q[0];
sx q[0];
rz(-0.73478886) q[0];
sx q[0];
rz(0.25811974) q[0];
rz(2.6737988) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(-1.3040257) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5060978) q[1];
sx q[1];
rz(-0.48759547) q[1];
sx q[1];
rz(2.0388842) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2673685) q[3];
sx q[3];
rz(-1.0099353) q[3];
sx q[3];
rz(-1.107533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.530431) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(-2.2116057) q[2];
rz(1.2215349) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4596443) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(-3.0514858) q[0];
rz(1.7168761) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(-0.43509126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6259436) q[0];
sx q[0];
rz(-2.0688147) q[0];
sx q[0];
rz(1.4048715) q[0];
x q[1];
rz(1.7516364) q[2];
sx q[2];
rz(-0.9107843) q[2];
sx q[2];
rz(1.7969799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96174091) q[1];
sx q[1];
rz(-1.6992178) q[1];
sx q[1];
rz(0.95203103) q[1];
rz(-pi) q[2];
rz(-1.5699205) q[3];
sx q[3];
rz(-1.9951576) q[3];
sx q[3];
rz(-0.76155886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-2.5153861) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(3.0873121) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(1.3351701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8568118) q[0];
sx q[0];
rz(-2.3216212) q[0];
sx q[0];
rz(-0.7818082) q[0];
x q[1];
rz(-0.35690709) q[2];
sx q[2];
rz(-1.0400598) q[2];
sx q[2];
rz(-0.70009795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21653342) q[1];
sx q[1];
rz(-2.0833942) q[1];
sx q[1];
rz(-0.61116265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.068633462) q[3];
sx q[3];
rz(-1.4331927) q[3];
sx q[3];
rz(-2.8559764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52428952) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-2.7900901) q[2];
rz(1.3737804) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-0.26941776) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3897301) q[0];
sx q[0];
rz(-0.26079145) q[0];
sx q[0];
rz(0.75916284) q[0];
rz(-2.0579386) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(-2.0738475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7368374) q[0];
sx q[0];
rz(-0.99946293) q[0];
sx q[0];
rz(2.570117) q[0];
x q[1];
rz(1.2933838) q[2];
sx q[2];
rz(-1.8545863) q[2];
sx q[2];
rz(-2.0744155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2764923) q[1];
sx q[1];
rz(-0.70521627) q[1];
sx q[1];
rz(0.92030763) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33371933) q[3];
sx q[3];
rz(-2.2401143) q[3];
sx q[3];
rz(-3.1143318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(-2.9313226) q[2];
rz(-2.353904) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44337153) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.5086077) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(-0.88268139) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(1.0630332) q[3];
sx q[3];
rz(-0.31217839) q[3];
sx q[3];
rz(0.97806539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
