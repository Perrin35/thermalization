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
rz(0.87585706) q[0];
sx q[0];
rz(-0.96867222) q[0];
sx q[0];
rz(-0.92585603) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449555) q[0];
sx q[0];
rz(-2.5384266) q[0];
sx q[0];
rz(-0.11267333) q[0];
x q[1];
rz(0.1466897) q[2];
sx q[2];
rz(-0.58977276) q[2];
sx q[2];
rz(1.1663933) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1086667) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(-1.0125748) q[1];
rz(-2.376287) q[3];
sx q[3];
rz(-1.0648784) q[3];
sx q[3];
rz(1.7047061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2609743) q[2];
sx q[2];
rz(-1.3857434) q[2];
sx q[2];
rz(-2.3495038) q[2];
rz(-0.875862) q[3];
sx q[3];
rz(-3.0328817) q[3];
sx q[3];
rz(-0.66332269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6302014) q[0];
sx q[0];
rz(-1.317861) q[0];
sx q[0];
rz(0.9414916) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-0.082854465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64697826) q[0];
sx q[0];
rz(-2.4002693) q[0];
sx q[0];
rz(-0.14880081) q[0];
rz(-2.156636) q[2];
sx q[2];
rz(-0.42030605) q[2];
sx q[2];
rz(0.46520761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8585457) q[1];
sx q[1];
rz(-2.0783922) q[1];
sx q[1];
rz(-2.0920442) q[1];
x q[2];
rz(-2.3979009) q[3];
sx q[3];
rz(-2.3451248) q[3];
sx q[3];
rz(1.9232779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2443709) q[2];
sx q[2];
rz(-1.7873849) q[2];
sx q[2];
rz(1.8841891) q[2];
rz(2.2952378) q[3];
sx q[3];
rz(-0.1592764) q[3];
sx q[3];
rz(-2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988679) q[0];
sx q[0];
rz(-1.6833479) q[0];
sx q[0];
rz(-0.22920907) q[0];
rz(0.21408679) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(-2.7600938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70951916) q[0];
sx q[0];
rz(-2.1901178) q[0];
sx q[0];
rz(1.1362057) q[0];
rz(-0.19241706) q[2];
sx q[2];
rz(-2.2936818) q[2];
sx q[2];
rz(-1.3012127) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22893045) q[1];
sx q[1];
rz(-0.060256392) q[1];
sx q[1];
rz(-2.069239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34708628) q[3];
sx q[3];
rz(-1.5659837) q[3];
sx q[3];
rz(-2.6322685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4182338) q[2];
sx q[2];
rz(-2.8853719) q[2];
sx q[2];
rz(0.56824938) q[2];
rz(-1.7226284) q[3];
sx q[3];
rz(-1.4293554) q[3];
sx q[3];
rz(1.0770146) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028246183) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(0.25076732) q[0];
rz(0.084550683) q[1];
sx q[1];
rz(-1.0944347) q[1];
sx q[1];
rz(1.2006522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47287175) q[0];
sx q[0];
rz(-2.3568912) q[0];
sx q[0];
rz(1.0502104) q[0];
rz(-pi) q[1];
rz(-1.0798825) q[2];
sx q[2];
rz(-0.49188313) q[2];
sx q[2];
rz(-1.5747923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29087092) q[1];
sx q[1];
rz(-0.71858908) q[1];
sx q[1];
rz(0.1512089) q[1];
rz(-pi) q[2];
rz(-0.26504604) q[3];
sx q[3];
rz(-1.0599905) q[3];
sx q[3];
rz(3.0731415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-2.0189221) q[2];
sx q[2];
rz(1.8505081) q[2];
rz(2.71991) q[3];
sx q[3];
rz(-2.6899874) q[3];
sx q[3];
rz(-1.5765367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6186433) q[0];
sx q[0];
rz(-2.4186501) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(0.067642637) q[1];
sx q[1];
rz(-1.5539955) q[1];
sx q[1];
rz(-1.1281475) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371409) q[0];
sx q[0];
rz(-1.6078976) q[0];
sx q[0];
rz(-0.18795975) q[0];
rz(-pi) q[1];
rz(1.0063926) q[2];
sx q[2];
rz(-2.5395406) q[2];
sx q[2];
rz(2.4133701) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50666821) q[1];
sx q[1];
rz(-1.4574058) q[1];
sx q[1];
rz(3.0874814) q[1];
rz(-pi) q[2];
rz(-0.87781436) q[3];
sx q[3];
rz(-2.7846309) q[3];
sx q[3];
rz(1.9138543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1077914) q[2];
sx q[2];
rz(-1.1148323) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(2.233861) q[3];
sx q[3];
rz(-1.9151442) q[3];
sx q[3];
rz(-0.41456732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229729) q[0];
sx q[0];
rz(-2.3899879) q[0];
sx q[0];
rz(-2.3722017) q[0];
rz(-1.9006624) q[1];
sx q[1];
rz(-0.85378328) q[1];
sx q[1];
rz(0.21533899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071393911) q[0];
sx q[0];
rz(-0.60766534) q[0];
sx q[0];
rz(1.2383226) q[0];
rz(0.017069503) q[2];
sx q[2];
rz(-1.9015549) q[2];
sx q[2];
rz(2.0278553) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78241888) q[1];
sx q[1];
rz(-2.4451588) q[1];
sx q[1];
rz(-0.38621186) q[1];
x q[2];
rz(1.2329383) q[3];
sx q[3];
rz(-2.2400899) q[3];
sx q[3];
rz(0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.874959) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(2.5249262) q[2];
rz(-0.2002317) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(-2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1277593) q[0];
sx q[0];
rz(-2.5795689) q[0];
sx q[0];
rz(3.1264937) q[0];
rz(-3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(2.1133568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5153582) q[0];
sx q[0];
rz(-1.6263824) q[0];
sx q[0];
rz(-0.65599156) q[0];
rz(-pi) q[1];
x q[1];
rz(2.188855) q[2];
sx q[2];
rz(-1.3985123) q[2];
sx q[2];
rz(-3.1069482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0575917) q[1];
sx q[1];
rz(-2.0023953) q[1];
sx q[1];
rz(-1.5398542) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3002197) q[3];
sx q[3];
rz(-1.8755873) q[3];
sx q[3];
rz(-2.6255325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5369109) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(-2.6417285) q[2];
rz(-0.31575051) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(1.0189112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71378088) q[0];
sx q[0];
rz(-1.2232057) q[0];
sx q[0];
rz(-0.94386238) q[0];
rz(1.7367412) q[1];
sx q[1];
rz(-1.2662788) q[1];
sx q[1];
rz(0.92165438) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0269146) q[0];
sx q[0];
rz(-1.252838) q[0];
sx q[0];
rz(1.5131895) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8878288) q[2];
sx q[2];
rz(-1.8602588) q[2];
sx q[2];
rz(-0.34115215) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4186395) q[1];
sx q[1];
rz(-2.2721842) q[1];
sx q[1];
rz(-2.4469814) q[1];
x q[2];
rz(1.5632024) q[3];
sx q[3];
rz(-0.51285997) q[3];
sx q[3];
rz(0.22554413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9965808) q[2];
sx q[2];
rz(-1.2715481) q[2];
sx q[2];
rz(1.5865631) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.8127706) q[3];
sx q[3];
rz(1.4012977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911875) q[0];
sx q[0];
rz(-2.0441971) q[0];
sx q[0];
rz(-0.64252585) q[0];
rz(1.4379028) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(-0.14370758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7449581) q[0];
sx q[0];
rz(-1.8734583) q[0];
sx q[0];
rz(0.58505262) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2238067) q[2];
sx q[2];
rz(-2.3812903) q[2];
sx q[2];
rz(1.1117293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4381126) q[1];
sx q[1];
rz(-0.92877711) q[1];
sx q[1];
rz(1.9025406) q[1];
x q[2];
rz(-1.1863352) q[3];
sx q[3];
rz(-1.1610306) q[3];
sx q[3];
rz(2.055197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6045427) q[2];
sx q[2];
rz(-1.1979251) q[2];
sx q[2];
rz(-0.26270467) q[2];
rz(-2.883834) q[3];
sx q[3];
rz(-2.7596605) q[3];
sx q[3];
rz(-2.2885382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8464552) q[0];
sx q[0];
rz(-1.4895804) q[0];
sx q[0];
rz(-1.2204131) q[0];
rz(2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(-0.89734546) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663956) q[0];
sx q[0];
rz(-2.0854073) q[0];
sx q[0];
rz(-1.4304016) q[0];
x q[1];
rz(0.32063516) q[2];
sx q[2];
rz(-2.8047987) q[2];
sx q[2];
rz(-1.8048665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7792977) q[1];
sx q[1];
rz(-0.26302281) q[1];
sx q[1];
rz(1.3785326) q[1];
rz(-pi) q[2];
rz(-1.4225716) q[3];
sx q[3];
rz(-2.5276042) q[3];
sx q[3];
rz(0.84166354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4974978) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(0.11478718) q[2];
rz(-0.74350205) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(-0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153618) q[0];
sx q[0];
rz(-0.76114934) q[0];
sx q[0];
rz(2.1834955) q[0];
rz(-1.505898) q[1];
sx q[1];
rz(-1.1264569) q[1];
sx q[1];
rz(1.1503848) q[1];
rz(-0.68810473) q[2];
sx q[2];
rz(-2.2672014) q[2];
sx q[2];
rz(-1.3939569) q[2];
rz(-0.43396797) q[3];
sx q[3];
rz(-1.4186191) q[3];
sx q[3];
rz(-0.71604244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
