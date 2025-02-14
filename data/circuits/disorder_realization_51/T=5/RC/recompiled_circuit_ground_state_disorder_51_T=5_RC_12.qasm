OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.12282523) q[0];
sx q[0];
rz(-0.0039761877) q[0];
sx q[0];
rz(-3.0206326) q[0];
rz(-2.6236985) q[1];
sx q[1];
rz(-2.8105812) q[1];
sx q[1];
rz(-1.9537227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5802104) q[0];
sx q[0];
rz(-0.31162173) q[0];
sx q[0];
rz(-2.8557957) q[0];
x q[1];
rz(-2.9953674) q[2];
sx q[2];
rz(-1.0492965) q[2];
sx q[2];
rz(-2.5086049) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49379865) q[1];
sx q[1];
rz(-2.0980853) q[1];
sx q[1];
rz(-2.6721558) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5224321) q[3];
sx q[3];
rz(-0.93975583) q[3];
sx q[3];
rz(-3.0121355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0594222) q[2];
sx q[2];
rz(-0.94258451) q[2];
sx q[2];
rz(1.7131294) q[2];
rz(-1.3905455) q[3];
sx q[3];
rz(-0.7775375) q[3];
sx q[3];
rz(-2.9318504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9222337) q[0];
sx q[0];
rz(-1.5809504) q[0];
sx q[0];
rz(-1.4571762) q[0];
rz(-0.68151418) q[1];
sx q[1];
rz(-2.4512955) q[1];
sx q[1];
rz(2.1843279) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0326704) q[0];
sx q[0];
rz(-1.7764979) q[0];
sx q[0];
rz(-2.6673596) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8419015) q[2];
sx q[2];
rz(-0.77163358) q[2];
sx q[2];
rz(-0.51811213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4034526) q[1];
sx q[1];
rz(-2.7305718) q[1];
sx q[1];
rz(-1.4477801) q[1];
rz(1.718347) q[3];
sx q[3];
rz(-3.0719764) q[3];
sx q[3];
rz(0.92629647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35931453) q[2];
sx q[2];
rz(-1.0474397) q[2];
sx q[2];
rz(-2.2347343) q[2];
rz(0.67306972) q[3];
sx q[3];
rz(-0.38452092) q[3];
sx q[3];
rz(1.8460021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68536711) q[0];
sx q[0];
rz(-2.4816368) q[0];
sx q[0];
rz(-0.54687706) q[0];
rz(1.3705672) q[1];
sx q[1];
rz(-1.1616881) q[1];
sx q[1];
rz(-3.0414157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4456383) q[0];
sx q[0];
rz(-1.4379408) q[0];
sx q[0];
rz(-1.8376875) q[0];
rz(-pi) q[1];
rz(-0.50649001) q[2];
sx q[2];
rz(-1.2236986) q[2];
sx q[2];
rz(0.60817761) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0292058) q[1];
sx q[1];
rz(-0.76923623) q[1];
sx q[1];
rz(2.8112414) q[1];
rz(-0.46296316) q[3];
sx q[3];
rz(-1.1826666) q[3];
sx q[3];
rz(2.8007647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3935516) q[2];
sx q[2];
rz(-1.2057722) q[2];
sx q[2];
rz(1.0432517) q[2];
rz(0.68759632) q[3];
sx q[3];
rz(-0.50191534) q[3];
sx q[3];
rz(-2.8858394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2954751) q[0];
sx q[0];
rz(-2.6961374) q[0];
sx q[0];
rz(1.0365781) q[0];
rz(-1.2713894) q[1];
sx q[1];
rz(-2.3174353) q[1];
sx q[1];
rz(1.2870671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4743253) q[0];
sx q[0];
rz(-1.0762574) q[0];
sx q[0];
rz(-0.70566341) q[0];
rz(-0.66133581) q[2];
sx q[2];
rz(-1.2244389) q[2];
sx q[2];
rz(-2.6544184) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85216552) q[1];
sx q[1];
rz(-2.1680538) q[1];
sx q[1];
rz(-2.7911623) q[1];
rz(-pi) q[2];
rz(-1.1274476) q[3];
sx q[3];
rz(-1.8487159) q[3];
sx q[3];
rz(-1.9390566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6048772) q[2];
sx q[2];
rz(-0.36467364) q[2];
sx q[2];
rz(3.1216915) q[2];
rz(-0.59602916) q[3];
sx q[3];
rz(-2.8607131) q[3];
sx q[3];
rz(1.178044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6562281) q[0];
sx q[0];
rz(-1.8940268) q[0];
sx q[0];
rz(-1.9670991) q[0];
rz(-2.8354722) q[1];
sx q[1];
rz(-1.0447634) q[1];
sx q[1];
rz(1.7721446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6868389) q[0];
sx q[0];
rz(-3.1343705) q[0];
sx q[0];
rz(2.0373086) q[0];
x q[1];
rz(-1.3197863) q[2];
sx q[2];
rz(-2.5282871) q[2];
sx q[2];
rz(1.5410739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8597653) q[1];
sx q[1];
rz(-1.2272115) q[1];
sx q[1];
rz(-0.1600034) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51403127) q[3];
sx q[3];
rz(-1.2054878) q[3];
sx q[3];
rz(-2.0333872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8451763) q[2];
sx q[2];
rz(-2.3974819) q[2];
sx q[2];
rz(-2.3386492) q[2];
rz(2.0261649) q[3];
sx q[3];
rz(-0.69474703) q[3];
sx q[3];
rz(-1.7723134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1882316) q[0];
sx q[0];
rz(-2.6022434) q[0];
sx q[0];
rz(-3.0292335) q[0];
rz(-2.864783) q[1];
sx q[1];
rz(-1.1667629) q[1];
sx q[1];
rz(-2.4940431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1153664) q[0];
sx q[0];
rz(-1.0293959) q[0];
sx q[0];
rz(2.7458835) q[0];
rz(2.0306088) q[2];
sx q[2];
rz(-2.1868878) q[2];
sx q[2];
rz(1.7700218) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4989483) q[1];
sx q[1];
rz(-2.0181839) q[1];
sx q[1];
rz(1.4397463) q[1];
x q[2];
rz(0.1754976) q[3];
sx q[3];
rz(-1.8030522) q[3];
sx q[3];
rz(1.0009223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67927805) q[2];
sx q[2];
rz(-2.9913112) q[2];
sx q[2];
rz(1.798299) q[2];
rz(-0.51370931) q[3];
sx q[3];
rz(-1.6535583) q[3];
sx q[3];
rz(2.546379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7996456) q[0];
sx q[0];
rz(-1.3222313) q[0];
sx q[0];
rz(-0.42022002) q[0];
rz(1.7814024) q[1];
sx q[1];
rz(-1.1667292) q[1];
sx q[1];
rz(-2.3374048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4046832) q[0];
sx q[0];
rz(-2.1054287) q[0];
sx q[0];
rz(-2.695309) q[0];
rz(-pi) q[1];
rz(2.8560195) q[2];
sx q[2];
rz(-2.336506) q[2];
sx q[2];
rz(-0.4090763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43257624) q[1];
sx q[1];
rz(-2.1423999) q[1];
sx q[1];
rz(-0.17523058) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10905091) q[3];
sx q[3];
rz(-1.7348467) q[3];
sx q[3];
rz(-1.4993639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5516025) q[2];
sx q[2];
rz(-1.0162105) q[2];
sx q[2];
rz(-1.972398) q[2];
rz(-2.9679838) q[3];
sx q[3];
rz(-0.9089402) q[3];
sx q[3];
rz(-1.9420697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6731897) q[0];
sx q[0];
rz(-2.115374) q[0];
sx q[0];
rz(1.2029368) q[0];
rz(0.24083336) q[1];
sx q[1];
rz(-0.92643654) q[1];
sx q[1];
rz(0.9333207) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77887452) q[0];
sx q[0];
rz(-1.9068471) q[0];
sx q[0];
rz(-1.1707992) q[0];
x q[1];
rz(-0.9983409) q[2];
sx q[2];
rz(-1.3476552) q[2];
sx q[2];
rz(-0.23912341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.62040239) q[1];
sx q[1];
rz(-2.7082084) q[1];
sx q[1];
rz(0.49232011) q[1];
x q[2];
rz(2.1347789) q[3];
sx q[3];
rz(-1.6996592) q[3];
sx q[3];
rz(-0.098524898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0241218) q[2];
sx q[2];
rz(-1.4535707) q[2];
sx q[2];
rz(0.13548279) q[2];
rz(0.99902117) q[3];
sx q[3];
rz(-2.8935367) q[3];
sx q[3];
rz(-0.55618709) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67021543) q[0];
sx q[0];
rz(-2.0129634) q[0];
sx q[0];
rz(1.3742597) q[0];
rz(-1.0820092) q[1];
sx q[1];
rz(-2.3550985) q[1];
sx q[1];
rz(-1.5622004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53331965) q[0];
sx q[0];
rz(-2.1100304) q[0];
sx q[0];
rz(3.1031552) q[0];
rz(-pi) q[1];
rz(-2.5815046) q[2];
sx q[2];
rz(-1.9011407) q[2];
sx q[2];
rz(-2.7719088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1154707) q[1];
sx q[1];
rz(-1.2392095) q[1];
sx q[1];
rz(-2.0381171) q[1];
x q[2];
rz(2.2213577) q[3];
sx q[3];
rz(-1.2830858) q[3];
sx q[3];
rz(-0.77863151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18206) q[2];
sx q[2];
rz(-1.5740732) q[2];
sx q[2];
rz(-2.5173397) q[2];
rz(2.4652081) q[3];
sx q[3];
rz(-1.1712149) q[3];
sx q[3];
rz(-2.306365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(-0.71037173) q[0];
sx q[0];
rz(-1.2483163) q[0];
sx q[0];
rz(-1.7927908) q[0];
rz(0.30385941) q[1];
sx q[1];
rz(-1.6108578) q[1];
sx q[1];
rz(2.2644728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6346587) q[0];
sx q[0];
rz(-1.4262232) q[0];
sx q[0];
rz(2.4094482) q[0];
rz(-pi) q[1];
rz(1.2459846) q[2];
sx q[2];
rz(-1.4673675) q[2];
sx q[2];
rz(1.6425486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7925229) q[1];
sx q[1];
rz(-1.190509) q[1];
sx q[1];
rz(2.040634) q[1];
x q[2];
rz(0.95939674) q[3];
sx q[3];
rz(-1.093321) q[3];
sx q[3];
rz(1.9408664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0861686) q[2];
sx q[2];
rz(-2.7358027) q[2];
sx q[2];
rz(-2.1093192) q[2];
rz(-2.0530733) q[3];
sx q[3];
rz(-1.4884596) q[3];
sx q[3];
rz(-2.3672339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.1879723) q[0];
sx q[0];
rz(-2.3379876) q[0];
sx q[0];
rz(0.048820989) q[0];
rz(1.8232518) q[1];
sx q[1];
rz(-1.7346458) q[1];
sx q[1];
rz(0.55191747) q[1];
rz(2.7248028) q[2];
sx q[2];
rz(-1.8183397) q[2];
sx q[2];
rz(0.37879735) q[2];
rz(-1.134675) q[3];
sx q[3];
rz(-1.5584617) q[3];
sx q[3];
rz(-0.8212318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
