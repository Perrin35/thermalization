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
rz(1.4242564) q[0];
sx q[0];
rz(-1.2196701) q[0];
sx q[0];
rz(-0.57504672) q[0];
rz(-2.9991034) q[1];
sx q[1];
rz(-1.2187076) q[1];
sx q[1];
rz(-2.277318) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7413717) q[0];
sx q[0];
rz(-1.3793886) q[0];
sx q[0];
rz(1.7975397) q[0];
x q[1];
rz(-2.9881485) q[2];
sx q[2];
rz(-0.94509238) q[2];
sx q[2];
rz(-0.99391711) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4886538) q[1];
sx q[1];
rz(-1.0407742) q[1];
sx q[1];
rz(-1.5731792) q[1];
rz(-2.0686758) q[3];
sx q[3];
rz(-1.0607294) q[3];
sx q[3];
rz(0.74397457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.064726) q[2];
sx q[2];
rz(-1.6419502) q[2];
sx q[2];
rz(1.0092674) q[2];
rz(-0.81389728) q[3];
sx q[3];
rz(-0.32865694) q[3];
sx q[3];
rz(2.8024659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547884) q[0];
sx q[0];
rz(-1.7330994) q[0];
sx q[0];
rz(-2.8993697) q[0];
rz(-1.9107266) q[1];
sx q[1];
rz(-0.63175285) q[1];
sx q[1];
rz(0.79808527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4674007) q[0];
sx q[0];
rz(-2.5134006) q[0];
sx q[0];
rz(-0.19285658) q[0];
x q[1];
rz(-1.961444) q[2];
sx q[2];
rz(-0.36875781) q[2];
sx q[2];
rz(-2.7386896) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.68951571) q[1];
sx q[1];
rz(-0.59714666) q[1];
sx q[1];
rz(-2.4680016) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7834218) q[3];
sx q[3];
rz(-1.7209952) q[3];
sx q[3];
rz(-1.1800486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4855087) q[2];
sx q[2];
rz(-1.069243) q[2];
sx q[2];
rz(-0.83267027) q[2];
rz(2.5101856) q[3];
sx q[3];
rz(-2.4000945) q[3];
sx q[3];
rz(-3.1413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176508) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(-0.18774524) q[0];
rz(2.2979157) q[1];
sx q[1];
rz(-1.2742821) q[1];
sx q[1];
rz(-0.63293308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62719856) q[0];
sx q[0];
rz(-2.7717675) q[0];
sx q[0];
rz(-1.7677977) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68667163) q[2];
sx q[2];
rz(-1.8455659) q[2];
sx q[2];
rz(-0.49897721) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5714076) q[1];
sx q[1];
rz(-1.8708036) q[1];
sx q[1];
rz(-0.5350929) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0121481) q[3];
sx q[3];
rz(-1.9151546) q[3];
sx q[3];
rz(0.5705006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26744276) q[2];
sx q[2];
rz(-2.1953857) q[2];
sx q[2];
rz(1.2713185) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(-0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8895759) q[0];
sx q[0];
rz(-0.75279623) q[0];
sx q[0];
rz(0.65943199) q[0];
rz(1.8239498) q[1];
sx q[1];
rz(-1.8023856) q[1];
sx q[1];
rz(-2.3514294) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20445828) q[0];
sx q[0];
rz(-0.35872981) q[0];
sx q[0];
rz(-0.35939868) q[0];
x q[1];
rz(1.6035346) q[2];
sx q[2];
rz(-1.590609) q[2];
sx q[2];
rz(1.1184426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.850833) q[1];
sx q[1];
rz(-1.5052126) q[1];
sx q[1];
rz(-2.6640434) q[1];
rz(-pi) q[2];
rz(-1.843256) q[3];
sx q[3];
rz(-1.0333697) q[3];
sx q[3];
rz(-0.3934653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(-2.8455632) q[2];
rz(2.5791903) q[3];
sx q[3];
rz(-2.2735169) q[3];
sx q[3];
rz(-0.11399046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89428467) q[0];
sx q[0];
rz(-1.9534651) q[0];
sx q[0];
rz(-2.9344015) q[0];
rz(-1.0317135) q[1];
sx q[1];
rz(-1.1528015) q[1];
sx q[1];
rz(-1.656104) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88865796) q[0];
sx q[0];
rz(-0.91080785) q[0];
sx q[0];
rz(2.071991) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8515622) q[2];
sx q[2];
rz(-1.0286912) q[2];
sx q[2];
rz(-2.5190767) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.42669) q[1];
sx q[1];
rz(-1.6985006) q[1];
sx q[1];
rz(-0.35121484) q[1];
rz(-0.10778963) q[3];
sx q[3];
rz(-2.2985085) q[3];
sx q[3];
rz(0.58290304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.038387211) q[2];
sx q[2];
rz(-2.3296671) q[2];
sx q[2];
rz(-2.0337598) q[2];
rz(-0.92249089) q[3];
sx q[3];
rz(-1.731571) q[3];
sx q[3];
rz(0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6237727) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(0.22860953) q[0];
rz(2.8672583) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(2.6913604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666445) q[0];
sx q[0];
rz(-0.83047509) q[0];
sx q[0];
rz(-1.9148097) q[0];
x q[1];
rz(-1.3415496) q[2];
sx q[2];
rz(-1.563213) q[2];
sx q[2];
rz(-2.6784865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5094368) q[1];
sx q[1];
rz(-0.56427279) q[1];
sx q[1];
rz(-2.6585101) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0110596) q[3];
sx q[3];
rz(-0.478906) q[3];
sx q[3];
rz(-2.4541209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7776514) q[2];
sx q[2];
rz(-2.5886017) q[2];
sx q[2];
rz(0.027776329) q[2];
rz(0.76644301) q[3];
sx q[3];
rz(-1.981363) q[3];
sx q[3];
rz(0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22208333) q[0];
sx q[0];
rz(-1.6351901) q[0];
sx q[0];
rz(0.62244225) q[0];
rz(0.70603236) q[1];
sx q[1];
rz(-0.89203867) q[1];
sx q[1];
rz(0.58696729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76510274) q[0];
sx q[0];
rz(-2.4060632) q[0];
sx q[0];
rz(-0.88659783) q[0];
x q[1];
rz(-2.6640011) q[2];
sx q[2];
rz(-2.6048663) q[2];
sx q[2];
rz(1.1004741) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.40499726) q[1];
sx q[1];
rz(-2.6301503) q[1];
sx q[1];
rz(-0.14528017) q[1];
rz(0.8587647) q[3];
sx q[3];
rz(-1.5786577) q[3];
sx q[3];
rz(2.238439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.075835) q[2];
sx q[2];
rz(-1.7087874) q[2];
sx q[2];
rz(2.5471121) q[2];
rz(-2.8822656) q[3];
sx q[3];
rz(-1.8766873) q[3];
sx q[3];
rz(-1.2172786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038430564) q[0];
sx q[0];
rz(-2.2362464) q[0];
sx q[0];
rz(-0.578798) q[0];
rz(1.831306) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(-0.8126747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90234251) q[0];
sx q[0];
rz(-2.6391811) q[0];
sx q[0];
rz(2.8397296) q[0];
x q[1];
rz(1.9695639) q[2];
sx q[2];
rz(-2.5426455) q[2];
sx q[2];
rz(-2.2676165) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1393237) q[1];
sx q[1];
rz(-1.2096757) q[1];
sx q[1];
rz(-2.4206672) q[1];
x q[2];
rz(2.1401494) q[3];
sx q[3];
rz(-2.0096931) q[3];
sx q[3];
rz(2.2341408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7909214) q[2];
sx q[2];
rz(-1.8427883) q[2];
sx q[2];
rz(-2.217963) q[2];
rz(-1.0203429) q[3];
sx q[3];
rz(-0.89710051) q[3];
sx q[3];
rz(3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807121) q[0];
sx q[0];
rz(-1.4360282) q[0];
sx q[0];
rz(1.6545779) q[0];
rz(-3.0912073) q[1];
sx q[1];
rz(-2.3245508) q[1];
sx q[1];
rz(-0.64687669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85031434) q[0];
sx q[0];
rz(-0.76843699) q[0];
sx q[0];
rz(0.9356229) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0912618) q[2];
sx q[2];
rz(-2.1796436) q[2];
sx q[2];
rz(-0.11907585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7868294) q[1];
sx q[1];
rz(-1.1656875) q[1];
sx q[1];
rz(-1.7988324) q[1];
rz(-1.7842125) q[3];
sx q[3];
rz(-1.4239256) q[3];
sx q[3];
rz(1.1606806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50463027) q[2];
sx q[2];
rz(-1.5865822) q[2];
sx q[2];
rz(-0.75616765) q[2];
rz(1.6715096) q[3];
sx q[3];
rz(-0.78592891) q[3];
sx q[3];
rz(2.3362931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037755448) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(1.1081498) q[0];
rz(0.26421079) q[1];
sx q[1];
rz(-1.4690396) q[1];
sx q[1];
rz(-2.5392551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717404) q[0];
sx q[0];
rz(-0.30769545) q[0];
sx q[0];
rz(1.0956531) q[0];
rz(-pi) q[1];
rz(0.17669038) q[2];
sx q[2];
rz(-1.7222705) q[2];
sx q[2];
rz(3.1164411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8089167) q[1];
sx q[1];
rz(-2.4530468) q[1];
sx q[1];
rz(1.6469514) q[1];
x q[2];
rz(0.29903166) q[3];
sx q[3];
rz(-0.67291884) q[3];
sx q[3];
rz(-2.3367019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5084874) q[2];
sx q[2];
rz(-2.5999531) q[2];
sx q[2];
rz(2.2701021) q[2];
rz(1.9320711) q[3];
sx q[3];
rz(-0.28035527) q[3];
sx q[3];
rz(-0.59396321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7021983) q[0];
sx q[0];
rz(-1.3855423) q[0];
sx q[0];
rz(2.6227797) q[0];
rz(0.58204542) q[1];
sx q[1];
rz(-2.0379635) q[1];
sx q[1];
rz(1.4720974) q[1];
rz(-1.9459361) q[2];
sx q[2];
rz(-0.87004253) q[2];
sx q[2];
rz(-0.042849356) q[2];
rz(2.9357435) q[3];
sx q[3];
rz(-2.4424853) q[3];
sx q[3];
rz(1.427099) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
