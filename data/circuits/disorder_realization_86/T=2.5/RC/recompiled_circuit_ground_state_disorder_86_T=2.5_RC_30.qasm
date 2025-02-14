OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0107083) q[0];
sx q[0];
rz(3.8151725) q[0];
sx q[0];
rz(7.1060915) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075283) q[0];
sx q[0];
rz(-0.82053608) q[0];
sx q[0];
rz(2.0506917) q[0];
rz(-pi) q[1];
rz(-2.1291586) q[2];
sx q[2];
rz(-1.867138) q[2];
sx q[2];
rz(1.7425962) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.339141) q[1];
sx q[1];
rz(-1.4467518) q[1];
sx q[1];
rz(1.5088827) q[1];
rz(-pi) q[2];
rz(-1.3884329) q[3];
sx q[3];
rz(-0.37624255) q[3];
sx q[3];
rz(-2.0227423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9660008) q[2];
sx q[2];
rz(-1.5660183) q[2];
sx q[2];
rz(1.9123745) q[2];
rz(1.843533) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(1.1840597) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9640279) q[0];
sx q[0];
rz(-2.8819045) q[0];
sx q[0];
rz(2.2143256) q[0];
rz(-2.941653) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(0.98958579) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67538136) q[0];
sx q[0];
rz(-2.7635734) q[0];
sx q[0];
rz(1.2751352) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43439718) q[2];
sx q[2];
rz(-1.6100307) q[2];
sx q[2];
rz(-2.5920282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.687508) q[1];
sx q[1];
rz(-1.5224147) q[1];
sx q[1];
rz(-1.6475186) q[1];
rz(-pi) q[2];
rz(-2.5603676) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(3.0581829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(2.7878063) q[2];
rz(0.95156041) q[3];
sx q[3];
rz(-0.74717251) q[3];
sx q[3];
rz(0.32683867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3554409) q[0];
sx q[0];
rz(-2.1936301) q[0];
sx q[0];
rz(-2.1308664) q[0];
rz(0.058723681) q[1];
sx q[1];
rz(-1.5208236) q[1];
sx q[1];
rz(-0.40930632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053913) q[0];
sx q[0];
rz(-2.276148) q[0];
sx q[0];
rz(-0.58569293) q[0];
rz(-1.204973) q[2];
sx q[2];
rz(-2.7049148) q[2];
sx q[2];
rz(-2.9819289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.066591) q[1];
sx q[1];
rz(-1.3745054) q[1];
sx q[1];
rz(2.58954) q[1];
x q[2];
rz(0.50893754) q[3];
sx q[3];
rz(-1.5526315) q[3];
sx q[3];
rz(0.048232676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0851486) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(2.9819581) q[2];
rz(-1.2498445) q[3];
sx q[3];
rz(-1.7227453) q[3];
sx q[3];
rz(-1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1551664) q[0];
sx q[0];
rz(-0.43743375) q[0];
sx q[0];
rz(-2.0470108) q[0];
rz(-1.1711586) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(1.0135244) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0243624) q[0];
sx q[0];
rz(-1.6147531) q[0];
sx q[0];
rz(1.5323601) q[0];
rz(-pi) q[1];
rz(-0.0036448467) q[2];
sx q[2];
rz(-1.0763775) q[2];
sx q[2];
rz(-1.0629176) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5202433) q[1];
sx q[1];
rz(-1.831437) q[1];
sx q[1];
rz(-1.5288749) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7988241) q[3];
sx q[3];
rz(-0.44525351) q[3];
sx q[3];
rz(-2.5955615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0052428) q[2];
sx q[2];
rz(-1.760027) q[2];
sx q[2];
rz(1.3137777) q[2];
rz(-2.2037196) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90792847) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(1.044957) q[0];
rz(-0.16432556) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(-2.9716861) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2192046) q[0];
sx q[0];
rz(-1.6419852) q[0];
sx q[0];
rz(-1.3842596) q[0];
rz(-pi) q[1];
rz(-0.81973524) q[2];
sx q[2];
rz(-2.6728874) q[2];
sx q[2];
rz(0.021989487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.610534) q[1];
sx q[1];
rz(-1.3278759) q[1];
sx q[1];
rz(-2.2353735) q[1];
rz(1.702013) q[3];
sx q[3];
rz(-2.2926712) q[3];
sx q[3];
rz(3.1027681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56875151) q[2];
sx q[2];
rz(-2.4390287) q[2];
sx q[2];
rz(-1.0315726) q[2];
rz(1.0860363) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.731855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3373435) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.401249) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-0.88637561) q[1];
sx q[1];
rz(1.8362129) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84112924) q[0];
sx q[0];
rz(-2.3337737) q[0];
sx q[0];
rz(0.43231583) q[0];
rz(3.1370509) q[2];
sx q[2];
rz(-2.0097187) q[2];
sx q[2];
rz(2.2786841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3363477) q[1];
sx q[1];
rz(-1.8191531) q[1];
sx q[1];
rz(-1.5389561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3280198) q[3];
sx q[3];
rz(-1.7194028) q[3];
sx q[3];
rz(-2.7000303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(-2.1567832) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(0.28520939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8816836) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(-0.024070865) q[1];
sx q[1];
rz(-1.6203974) q[1];
sx q[1];
rz(-0.37758652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1771072) q[0];
sx q[0];
rz(-1.7269352) q[0];
sx q[0];
rz(0.6556169) q[0];
rz(-1.7837011) q[2];
sx q[2];
rz(-1.3981552) q[2];
sx q[2];
rz(1.6253234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3096034) q[1];
sx q[1];
rz(-1.1909232) q[1];
sx q[1];
rz(1.1516476) q[1];
x q[2];
rz(2.0041647) q[3];
sx q[3];
rz(-1.876017) q[3];
sx q[3];
rz(2.5030672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5923803) q[2];
sx q[2];
rz(-0.33752957) q[2];
sx q[2];
rz(2.737992) q[2];
rz(-0.18925439) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(-2.6846867) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6117578) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-0.31103617) q[0];
rz(2.1850736) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(-2.8172353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2431902) q[0];
sx q[0];
rz(-2.0721738) q[0];
sx q[0];
rz(1.7883975) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10559428) q[2];
sx q[2];
rz(-2.6289231) q[2];
sx q[2];
rz(1.6216506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2203103) q[1];
sx q[1];
rz(-1.9844354) q[1];
sx q[1];
rz(0.92343729) q[1];
x q[2];
rz(1.1127541) q[3];
sx q[3];
rz(-2.4107217) q[3];
sx q[3];
rz(0.79938408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4057464) q[2];
sx q[2];
rz(-2.3509071) q[2];
sx q[2];
rz(-0.22845593) q[2];
rz(0.53260803) q[3];
sx q[3];
rz(-0.76627982) q[3];
sx q[3];
rz(-2.2920091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5936977) q[0];
sx q[0];
rz(-2.6314647) q[0];
sx q[0];
rz(-1.8713895) q[0];
rz(-0.58865976) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(-2.7587845) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24341881) q[0];
sx q[0];
rz(-0.23056689) q[0];
sx q[0];
rz(-1.4617993) q[0];
rz(-1.8201939) q[2];
sx q[2];
rz(-1.9928842) q[2];
sx q[2];
rz(0.37427926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11483604) q[1];
sx q[1];
rz(-0.72163218) q[1];
sx q[1];
rz(0.19087692) q[1];
rz(-2.8715897) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(1.6635513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9855014) q[2];
sx q[2];
rz(-1.5727377) q[2];
sx q[2];
rz(-0.4001948) q[2];
rz(2.8943446) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8409214) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(2.4726942) q[1];
sx q[1];
rz(-1.4543507) q[1];
sx q[1];
rz(1.5060172) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6164613) q[0];
sx q[0];
rz(-2.2716953) q[0];
sx q[0];
rz(0.67411311) q[0];
rz(-2.593562) q[2];
sx q[2];
rz(-1.1924253) q[2];
sx q[2];
rz(-3.1032094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4989711) q[1];
sx q[1];
rz(-1.5318512) q[1];
sx q[1];
rz(-2.0508815) q[1];
x q[2];
rz(2.5745939) q[3];
sx q[3];
rz(-0.97494379) q[3];
sx q[3];
rz(-1.8351549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(1.6157185) q[2];
rz(0.29159355) q[3];
sx q[3];
rz(-2.6988131) q[3];
sx q[3];
rz(-0.35167882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780554) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(-0.28868227) q[1];
sx q[1];
rz(-2.3434227) q[1];
sx q[1];
rz(-3.1176288) q[1];
rz(2.504442) q[2];
sx q[2];
rz(-0.49497866) q[2];
sx q[2];
rz(-0.53866932) q[2];
rz(-0.45997672) q[3];
sx q[3];
rz(-1.3446076) q[3];
sx q[3];
rz(1.4190058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
