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
rz(0.54230827) q[0];
sx q[0];
rz(-0.13442726) q[0];
sx q[0];
rz(-1.0472714) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(2.9919762) q[1];
sx q[1];
rz(11.621578) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0243264) q[0];
sx q[0];
rz(-0.42792441) q[0];
sx q[0];
rz(-0.86366349) q[0];
rz(2.216835) q[2];
sx q[2];
rz(-1.7584718) q[2];
sx q[2];
rz(0.061522324) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8103247) q[1];
sx q[1];
rz(-2.4400418) q[1];
sx q[1];
rz(-2.3980106) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0534001) q[3];
sx q[3];
rz(-1.6783829) q[3];
sx q[3];
rz(1.3914598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2904539) q[2];
sx q[2];
rz(-1.4522499) q[2];
sx q[2];
rz(-1.5645082) q[2];
rz(2.5751298) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(1.2982781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-0.7199921) q[0];
rz(2.8672245) q[1];
sx q[1];
rz(-2.4101078) q[1];
sx q[1];
rz(-0.55363399) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1328515) q[0];
sx q[0];
rz(-2.7692502) q[0];
sx q[0];
rz(-1.8944064) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4963412) q[2];
sx q[2];
rz(-0.49075365) q[2];
sx q[2];
rz(0.09749271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95018847) q[1];
sx q[1];
rz(-0.85655115) q[1];
sx q[1];
rz(-0.49226239) q[1];
x q[2];
rz(-2.6194044) q[3];
sx q[3];
rz(-2.0117674) q[3];
sx q[3];
rz(-0.14312927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1729892) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(2.7776862) q[2];
rz(-0.11161741) q[3];
sx q[3];
rz(-2.2026187) q[3];
sx q[3];
rz(0.76044559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97742057) q[0];
sx q[0];
rz(-0.82472473) q[0];
sx q[0];
rz(0.0053996276) q[0];
rz(-2.3258356) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(0.38591787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49026981) q[0];
sx q[0];
rz(-1.2322591) q[0];
sx q[0];
rz(-1.8671892) q[0];
x q[1];
rz(-1.0089985) q[2];
sx q[2];
rz(-1.2616065) q[2];
sx q[2];
rz(-0.55436347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2009243) q[1];
sx q[1];
rz(-1.4438119) q[1];
sx q[1];
rz(1.9534355) q[1];
x q[2];
rz(1.98514) q[3];
sx q[3];
rz(-1.6938391) q[3];
sx q[3];
rz(-2.4407354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5526814) q[2];
sx q[2];
rz(-0.83325714) q[2];
sx q[2];
rz(-1.3820648) q[2];
rz(-2.9803993) q[3];
sx q[3];
rz(-1.7101945) q[3];
sx q[3];
rz(2.5126357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88518077) q[0];
sx q[0];
rz(-1.8653402) q[0];
sx q[0];
rz(-1.3389583) q[0];
rz(1.3171875) q[1];
sx q[1];
rz(-0.97511292) q[1];
sx q[1];
rz(-0.29979527) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3835589) q[0];
sx q[0];
rz(-1.9179317) q[0];
sx q[0];
rz(-0.68521037) q[0];
rz(-2.6671403) q[2];
sx q[2];
rz(-2.1079194) q[2];
sx q[2];
rz(-0.65309292) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2721709) q[1];
sx q[1];
rz(-1.4336839) q[1];
sx q[1];
rz(-1.9253982) q[1];
rz(-pi) q[2];
rz(2.2065877) q[3];
sx q[3];
rz(-1.4640558) q[3];
sx q[3];
rz(0.028117953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6224711) q[2];
sx q[2];
rz(-0.93834472) q[2];
sx q[2];
rz(2.8224714) q[2];
rz(3.0506813) q[3];
sx q[3];
rz(-0.14486434) q[3];
sx q[3];
rz(-1.7849785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4153445) q[0];
sx q[0];
rz(-0.81979668) q[0];
sx q[0];
rz(-3.1392642) q[0];
rz(1.5709411) q[1];
sx q[1];
rz(-0.80051533) q[1];
sx q[1];
rz(1.9470107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1462458) q[0];
sx q[0];
rz(-1.6888535) q[0];
sx q[0];
rz(0.15583584) q[0];
x q[1];
rz(1.5244687) q[2];
sx q[2];
rz(-0.71832685) q[2];
sx q[2];
rz(-0.83889942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71517005) q[1];
sx q[1];
rz(-1.914901) q[1];
sx q[1];
rz(3.032883) q[1];
rz(-1.8525328) q[3];
sx q[3];
rz(-2.1093371) q[3];
sx q[3];
rz(-1.0823712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.883256) q[2];
sx q[2];
rz(-0.77815762) q[2];
sx q[2];
rz(3.0641595) q[2];
rz(2.1954913) q[3];
sx q[3];
rz(-2.6587722) q[3];
sx q[3];
rz(-0.4536804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9282064) q[0];
sx q[0];
rz(-3.0551857) q[0];
sx q[0];
rz(1.7267831) q[0];
rz(0.66697031) q[1];
sx q[1];
rz(-1.7844776) q[1];
sx q[1];
rz(0.40474969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20807438) q[0];
sx q[0];
rz(-0.41941038) q[0];
sx q[0];
rz(2.389877) q[0];
rz(-pi) q[1];
rz(0.26164345) q[2];
sx q[2];
rz(-2.2115797) q[2];
sx q[2];
rz(1.4705758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0424543) q[1];
sx q[1];
rz(-2.6746866) q[1];
sx q[1];
rz(-1.7068638) q[1];
x q[2];
rz(2.1920106) q[3];
sx q[3];
rz(-0.57692617) q[3];
sx q[3];
rz(-0.85437894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82379597) q[2];
sx q[2];
rz(-2.3522289) q[2];
sx q[2];
rz(-2.6981567) q[2];
rz(2.7277842) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(-2.2558291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39171788) q[0];
sx q[0];
rz(-1.7789919) q[0];
sx q[0];
rz(2.474127) q[0];
rz(-0.04774566) q[1];
sx q[1];
rz(-1.1404488) q[1];
sx q[1];
rz(1.1183636) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76343173) q[0];
sx q[0];
rz(-1.1677647) q[0];
sx q[0];
rz(0.58149882) q[0];
rz(-1.7553545) q[2];
sx q[2];
rz(-0.58635752) q[2];
sx q[2];
rz(1.9826012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0065919) q[1];
sx q[1];
rz(-1.9891095) q[1];
sx q[1];
rz(0.63726421) q[1];
rz(1.4468071) q[3];
sx q[3];
rz(-0.50278864) q[3];
sx q[3];
rz(-3.0320252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.095801) q[2];
sx q[2];
rz(-1.8818776) q[2];
sx q[2];
rz(0.0087139159) q[2];
rz(1.2039315) q[3];
sx q[3];
rz(-0.34398505) q[3];
sx q[3];
rz(-3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0225723) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(-0.88985306) q[0];
rz(1.8564557) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(3.0056675) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0014627731) q[0];
sx q[0];
rz(-2.9076932) q[0];
sx q[0];
rz(2.9359096) q[0];
rz(2.1077584) q[2];
sx q[2];
rz(-1.0091678) q[2];
sx q[2];
rz(-1.461535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0697729) q[1];
sx q[1];
rz(-1.7679224) q[1];
sx q[1];
rz(-0.33918753) q[1];
x q[2];
rz(2.9985524) q[3];
sx q[3];
rz(-1.2238127) q[3];
sx q[3];
rz(0.59511371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70302427) q[2];
sx q[2];
rz(-0.88006222) q[2];
sx q[2];
rz(-2.0795889) q[2];
rz(0.92987531) q[3];
sx q[3];
rz(-0.78571856) q[3];
sx q[3];
rz(-0.78833956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.7787665) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(0.40089259) q[0];
rz(-0.76155424) q[1];
sx q[1];
rz(-1.4925894) q[1];
sx q[1];
rz(0.78899312) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74247201) q[0];
sx q[0];
rz(-1.443304) q[0];
sx q[0];
rz(1.5546868) q[0];
rz(-pi) q[1];
rz(-1.3775741) q[2];
sx q[2];
rz(-0.68422645) q[2];
sx q[2];
rz(-2.2598337) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4915062) q[1];
sx q[1];
rz(-1.6900926) q[1];
sx q[1];
rz(2.5507798) q[1];
rz(-pi) q[2];
rz(-0.085461334) q[3];
sx q[3];
rz(-2.6151867) q[3];
sx q[3];
rz(1.5554903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5679428) q[2];
sx q[2];
rz(-1.8486134) q[2];
sx q[2];
rz(2.4207777) q[2];
rz(1.6503664) q[3];
sx q[3];
rz(-0.78250116) q[3];
sx q[3];
rz(-1.8318374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023130527) q[0];
sx q[0];
rz(-3.0266302) q[0];
sx q[0];
rz(2.3638828) q[0];
rz(-2.481781) q[1];
sx q[1];
rz(-2.2751364) q[1];
sx q[1];
rz(0.39001098) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5725419) q[0];
sx q[0];
rz(-3.0016698) q[0];
sx q[0];
rz(-1.8403017) q[0];
rz(-pi) q[1];
rz(-0.74441461) q[2];
sx q[2];
rz(-1.0420024) q[2];
sx q[2];
rz(-2.2802558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9007787) q[1];
sx q[1];
rz(-1.1692746) q[1];
sx q[1];
rz(2.8177849) q[1];
rz(2.2220192) q[3];
sx q[3];
rz(-1.9863141) q[3];
sx q[3];
rz(-0.093747698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4593792) q[2];
sx q[2];
rz(-2.7028658) q[2];
sx q[2];
rz(-0.48975804) q[2];
rz(0.30719906) q[3];
sx q[3];
rz(-2.2178853) q[3];
sx q[3];
rz(-1.3908305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78140344) q[0];
sx q[0];
rz(-1.4908981) q[0];
sx q[0];
rz(1.4733462) q[0];
rz(-1.0745984) q[1];
sx q[1];
rz(-0.85292024) q[1];
sx q[1];
rz(-1.9180752) q[1];
rz(-2.3831153) q[2];
sx q[2];
rz(-0.52634326) q[2];
sx q[2];
rz(1.1545622) q[2];
rz(0.60081595) q[3];
sx q[3];
rz(-1.0944585) q[3];
sx q[3];
rz(-0.16437358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
