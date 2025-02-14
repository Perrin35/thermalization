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
rz(0.55573207) q[0];
sx q[0];
rz(4.421173) q[0];
sx q[0];
rz(9.0968994) q[0];
rz(0.15283395) q[1];
sx q[1];
rz(5.7938303) q[1];
sx q[1];
rz(7.2942776) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0450789) q[0];
sx q[0];
rz(-2.156684) q[0];
sx q[0];
rz(1.7631329) q[0];
x q[1];
rz(-1.7494124) q[2];
sx q[2];
rz(-1.9733323) q[2];
sx q[2];
rz(0.91031633) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0786405) q[1];
sx q[1];
rz(-0.83940369) q[1];
sx q[1];
rz(2.6432829) q[1];
rz(2.070363) q[3];
sx q[3];
rz(-0.46854436) q[3];
sx q[3];
rz(2.6234181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27935394) q[2];
sx q[2];
rz(-2.3591159) q[2];
sx q[2];
rz(1.5581101) q[2];
rz(0.33485788) q[3];
sx q[3];
rz(-2.0575276) q[3];
sx q[3];
rz(2.8607232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25664499) q[0];
sx q[0];
rz(-2.5769951) q[0];
sx q[0];
rz(-0.33194342) q[0];
rz(0.360082) q[1];
sx q[1];
rz(-1.8372476) q[1];
sx q[1];
rz(-0.28775451) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1126039) q[0];
sx q[0];
rz(-1.3224012) q[0];
sx q[0];
rz(1.364014) q[0];
rz(2.4096978) q[2];
sx q[2];
rz(-1.8127155) q[2];
sx q[2];
rz(-3.0223924) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.33720484) q[1];
sx q[1];
rz(-1.241695) q[1];
sx q[1];
rz(-1.2044524) q[1];
rz(-2.5310181) q[3];
sx q[3];
rz(-2.861851) q[3];
sx q[3];
rz(0.66204643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.687279) q[2];
sx q[2];
rz(-1.6087029) q[2];
sx q[2];
rz(-1.3909371) q[2];
rz(-2.2823997) q[3];
sx q[3];
rz(-1.2733368) q[3];
sx q[3];
rz(2.9505762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.745382) q[0];
sx q[0];
rz(-1.6062382) q[0];
sx q[0];
rz(-2.9111653) q[0];
rz(0.51741171) q[1];
sx q[1];
rz(-2.0473862) q[1];
sx q[1];
rz(2.8527625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.13076) q[0];
sx q[0];
rz(-1.9790181) q[0];
sx q[0];
rz(-0.064716332) q[0];
rz(-pi) q[1];
rz(-0.14634653) q[2];
sx q[2];
rz(-0.47292559) q[2];
sx q[2];
rz(2.2951916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9625712) q[1];
sx q[1];
rz(-1.5177392) q[1];
sx q[1];
rz(1.1607443) q[1];
rz(-pi) q[2];
rz(-1.3473116) q[3];
sx q[3];
rz(-1.399125) q[3];
sx q[3];
rz(0.61644256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1032224) q[2];
sx q[2];
rz(-2.4028845) q[2];
sx q[2];
rz(-0.040741097) q[2];
rz(-0.1013969) q[3];
sx q[3];
rz(-1.8056168) q[3];
sx q[3];
rz(-2.0749157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3876225) q[0];
sx q[0];
rz(-3.1356223) q[0];
sx q[0];
rz(2.907584) q[0];
rz(-0.19800828) q[1];
sx q[1];
rz(-2.104069) q[1];
sx q[1];
rz(-0.98349804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3458684) q[0];
sx q[0];
rz(-1.0875174) q[0];
sx q[0];
rz(-3.0935174) q[0];
x q[1];
rz(1.4794691) q[2];
sx q[2];
rz(-1.0589561) q[2];
sx q[2];
rz(0.78081607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3819921) q[1];
sx q[1];
rz(-0.97642094) q[1];
sx q[1];
rz(-0.71099736) q[1];
x q[2];
rz(1.4928905) q[3];
sx q[3];
rz(-0.39038218) q[3];
sx q[3];
rz(-1.4161033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6733751) q[2];
sx q[2];
rz(-0.28926352) q[2];
sx q[2];
rz(-0.77486983) q[2];
rz(-1.8153927) q[3];
sx q[3];
rz(-1.1893136) q[3];
sx q[3];
rz(-0.51138043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4020017) q[0];
sx q[0];
rz(-0.87257659) q[0];
sx q[0];
rz(-0.032489754) q[0];
rz(1.912311) q[1];
sx q[1];
rz(-0.928343) q[1];
sx q[1];
rz(-3.0585739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0262526) q[0];
sx q[0];
rz(-1.7949545) q[0];
sx q[0];
rz(-2.6782533) q[0];
rz(-pi) q[1];
x q[1];
rz(2.584721) q[2];
sx q[2];
rz(-2.8285714) q[2];
sx q[2];
rz(2.0890631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0159576) q[1];
sx q[1];
rz(-2.2803366) q[1];
sx q[1];
rz(2.5112391) q[1];
rz(-0.015542726) q[3];
sx q[3];
rz(-2.2076026) q[3];
sx q[3];
rz(-1.4525082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3518389) q[2];
sx q[2];
rz(-2.3408076) q[2];
sx q[2];
rz(2.9808673) q[2];
rz(-1.9488581) q[3];
sx q[3];
rz(-0.21218097) q[3];
sx q[3];
rz(0.76782697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47650325) q[0];
sx q[0];
rz(-1.1192717) q[0];
sx q[0];
rz(2.8506668) q[0];
rz(1.3833969) q[1];
sx q[1];
rz(-1.29888) q[1];
sx q[1];
rz(-0.49984041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4668149) q[0];
sx q[0];
rz(-1.5272015) q[0];
sx q[0];
rz(1.528426) q[0];
rz(-pi) q[1];
rz(2.2989474) q[2];
sx q[2];
rz(-1.1367186) q[2];
sx q[2];
rz(-1.7378716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0575057) q[1];
sx q[1];
rz(-1.1228787) q[1];
sx q[1];
rz(-2.491076) q[1];
x q[2];
rz(-0.53414102) q[3];
sx q[3];
rz(-2.2075966) q[3];
sx q[3];
rz(2.4780688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2073888) q[2];
sx q[2];
rz(-0.61177212) q[2];
sx q[2];
rz(-0.70154166) q[2];
rz(0.97405854) q[3];
sx q[3];
rz(-2.3249224) q[3];
sx q[3];
rz(0.90132236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262254) q[0];
sx q[0];
rz(-0.81804818) q[0];
sx q[0];
rz(2.4819964) q[0];
rz(-1.2391799) q[1];
sx q[1];
rz(-2.0678803) q[1];
sx q[1];
rz(-2.0674131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2240018) q[0];
sx q[0];
rz(-1.5573643) q[0];
sx q[0];
rz(-1.55634) q[0];
rz(-0.72788179) q[2];
sx q[2];
rz(-2.5932576) q[2];
sx q[2];
rz(2.9408583) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32666001) q[1];
sx q[1];
rz(-0.49100093) q[1];
sx q[1];
rz(-2.0424974) q[1];
rz(-pi) q[2];
rz(-0.76403615) q[3];
sx q[3];
rz(-1.8978528) q[3];
sx q[3];
rz(0.59248176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85019511) q[2];
sx q[2];
rz(-2.3829298) q[2];
sx q[2];
rz(-2.5568753) q[2];
rz(2.0753453) q[3];
sx q[3];
rz(-1.2277579) q[3];
sx q[3];
rz(0.40759459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9503815) q[0];
sx q[0];
rz(-1.1703015) q[0];
sx q[0];
rz(-0.48072746) q[0];
rz(1.2113384) q[1];
sx q[1];
rz(-1.6119266) q[1];
sx q[1];
rz(-0.617625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4803333) q[0];
sx q[0];
rz(-0.32513371) q[0];
sx q[0];
rz(-1.5839769) q[0];
x q[1];
rz(2.0938056) q[2];
sx q[2];
rz(-1.0289425) q[2];
sx q[2];
rz(-0.43719342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94719515) q[1];
sx q[1];
rz(-0.36809599) q[1];
sx q[1];
rz(0.93666623) q[1];
rz(1.6061574) q[3];
sx q[3];
rz(-1.3786982) q[3];
sx q[3];
rz(-2.0336322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9628613) q[2];
sx q[2];
rz(-1.3946673) q[2];
sx q[2];
rz(-2.2507131) q[2];
rz(1.0293055) q[3];
sx q[3];
rz(-1.3903214) q[3];
sx q[3];
rz(1.5911969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7903098) q[0];
sx q[0];
rz(-2.2064378) q[0];
sx q[0];
rz(-1.1736322) q[0];
rz(-1.952518) q[1];
sx q[1];
rz(-0.74780858) q[1];
sx q[1];
rz(0.48114166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3404589) q[0];
sx q[0];
rz(-2.1251571) q[0];
sx q[0];
rz(1.8340684) q[0];
rz(-pi) q[1];
rz(-1.6716154) q[2];
sx q[2];
rz(-1.013213) q[2];
sx q[2];
rz(2.4876311) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0125825) q[1];
sx q[1];
rz(-0.6193739) q[1];
sx q[1];
rz(-1.8257837) q[1];
rz(-pi) q[2];
rz(-2.2523426) q[3];
sx q[3];
rz(-0.64683952) q[3];
sx q[3];
rz(1.7689887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9159307) q[2];
sx q[2];
rz(-1.91232) q[2];
sx q[2];
rz(1.4813598) q[2];
rz(-3.0952752) q[3];
sx q[3];
rz(-1.9527083) q[3];
sx q[3];
rz(1.9143298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34761053) q[0];
sx q[0];
rz(-2.1078258) q[0];
sx q[0];
rz(0.069742918) q[0];
rz(2.8522885) q[1];
sx q[1];
rz(-1.8569088) q[1];
sx q[1];
rz(-2.0893673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5525411) q[0];
sx q[0];
rz(-0.49363401) q[0];
sx q[0];
rz(2.4049253) q[0];
rz(-pi) q[1];
rz(-2.9751189) q[2];
sx q[2];
rz(-2.2262339) q[2];
sx q[2];
rz(-0.85879281) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47202808) q[1];
sx q[1];
rz(-1.2109562) q[1];
sx q[1];
rz(1.1675952) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77731384) q[3];
sx q[3];
rz(-2.5717426) q[3];
sx q[3];
rz(0.22138813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59794402) q[2];
sx q[2];
rz(-1.9660549) q[2];
sx q[2];
rz(-3.0808466) q[2];
rz(-2.8525823) q[3];
sx q[3];
rz(-1.776639) q[3];
sx q[3];
rz(-1.2750767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.65171) q[0];
sx q[0];
rz(-1.4986421) q[0];
sx q[0];
rz(-1.0810252) q[0];
rz(1.0501077) q[1];
sx q[1];
rz(-1.0625912) q[1];
sx q[1];
rz(-1.2407632) q[1];
rz(-2.0403258) q[2];
sx q[2];
rz(-1.4360089) q[2];
sx q[2];
rz(0.77965005) q[2];
rz(-0.28239863) q[3];
sx q[3];
rz(-1.1544607) q[3];
sx q[3];
rz(-0.37731597) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
