OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(-2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(2.8957638) q[1];
sx q[1];
rz(10.413269) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91627097) q[0];
sx q[0];
rz(-2.7164408) q[0];
sx q[0];
rz(-1.7122373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97889401) q[2];
sx q[2];
rz(-1.4375293) q[2];
sx q[2];
rz(0.29083458) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.72585427) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(1.0456677) q[1];
rz(0.037976102) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(-1.2167041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6926379) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(2.5732178) q[2];
rz(-1.9521936) q[3];
sx q[3];
rz(-2.9247734) q[3];
sx q[3];
rz(2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73873591) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(0.96827132) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(1.2526858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90447146) q[0];
sx q[0];
rz(-1.4112817) q[0];
sx q[0];
rz(1.6617387) q[0];
x q[1];
rz(2.2717936) q[2];
sx q[2];
rz(-0.79248488) q[2];
sx q[2];
rz(-3.062641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0879841) q[1];
sx q[1];
rz(-1.0454185) q[1];
sx q[1];
rz(-2.5768075) q[1];
x q[2];
rz(-1.9545752) q[3];
sx q[3];
rz(-2.6740101) q[3];
sx q[3];
rz(0.13320696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12053717) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(-2.7978314) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-2.7407586) q[3];
sx q[3];
rz(0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47135982) q[0];
sx q[0];
rz(-0.25046644) q[0];
sx q[0];
rz(0.0070455889) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(0.71281707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757652) q[0];
sx q[0];
rz(-1.3607425) q[0];
sx q[0];
rz(-1.7957627) q[0];
x q[1];
rz(-2.9696434) q[2];
sx q[2];
rz(-1.9438582) q[2];
sx q[2];
rz(0.821515) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30046001) q[1];
sx q[1];
rz(-1.7736048) q[1];
sx q[1];
rz(-1.0293142) q[1];
rz(-pi) q[2];
rz(2.0531822) q[3];
sx q[3];
rz(-1.4995121) q[3];
sx q[3];
rz(-2.7977668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(2.4776069) q[2];
rz(1.9021696) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5666714) q[0];
sx q[0];
rz(-0.029733812) q[0];
sx q[0];
rz(-2.5675039) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(2.2132197) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61558047) q[0];
sx q[0];
rz(-0.86475879) q[0];
sx q[0];
rz(3.1403149) q[0];
rz(-pi) q[1];
rz(0.78903918) q[2];
sx q[2];
rz(-1.3256729) q[2];
sx q[2];
rz(-0.48691985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91765431) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(1.1397821) q[1];
rz(-0.5023605) q[3];
sx q[3];
rz(-2.1661048) q[3];
sx q[3];
rz(-2.3140698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.1553923) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-0.5345878) q[3];
sx q[3];
rz(2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(1.4986562) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(-2.1898851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18626801) q[0];
sx q[0];
rz(-0.82337472) q[0];
sx q[0];
rz(-3.0427409) q[0];
rz(-2.5809418) q[2];
sx q[2];
rz(-0.89386212) q[2];
sx q[2];
rz(1.4844984) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4068027) q[1];
sx q[1];
rz(-0.64195913) q[1];
sx q[1];
rz(-1.2924679) q[1];
rz(2.1390901) q[3];
sx q[3];
rz(-1.2270524) q[3];
sx q[3];
rz(-2.7078201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.1523694) q[2];
rz(2.0955829) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9933269) q[0];
sx q[0];
rz(-2.1874805) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(2.1173677) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(-1.5354059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.673296) q[0];
sx q[0];
rz(-1.6185986) q[0];
sx q[0];
rz(-3.0424546) q[0];
rz(-2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(-2.6256109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.225243) q[1];
sx q[1];
rz(-2.8147329) q[1];
sx q[1];
rz(-2.0845695) q[1];
rz(-pi) q[2];
rz(-0.225004) q[3];
sx q[3];
rz(-2.5000754) q[3];
sx q[3];
rz(-1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17807047) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(0.17573389) q[2];
rz(2.0641616) q[3];
sx q[3];
rz(-1.9947937) q[3];
sx q[3];
rz(-0.0065461672) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734633) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(1.8716795) q[0];
rz(0.1779671) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(1.7756745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1195927) q[0];
sx q[0];
rz(-2.3455142) q[0];
sx q[0];
rz(-0.47795602) q[0];
x q[1];
rz(-1.3283417) q[2];
sx q[2];
rz(-1.889466) q[2];
sx q[2];
rz(-2.6048425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3733702) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(-2.8819487) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9317022) q[3];
sx q[3];
rz(-2.6025747) q[3];
sx q[3];
rz(-0.15912661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0697249) q[2];
sx q[2];
rz(-0.29399997) q[2];
sx q[2];
rz(0.89795566) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(0.32178497) q[0];
rz(-0.92542648) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2935534) q[0];
sx q[0];
rz(-0.32395054) q[0];
sx q[0];
rz(2.241961) q[0];
x q[1];
rz(3.1084193) q[2];
sx q[2];
rz(-1.355194) q[2];
sx q[2];
rz(-2.6278091) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1083793) q[1];
sx q[1];
rz(-1.2011659) q[1];
sx q[1];
rz(2.0481678) q[1];
rz(1.6156322) q[3];
sx q[3];
rz(-2.3196844) q[3];
sx q[3];
rz(-2.3007948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-2.896538) q[3];
sx q[3];
rz(1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(1.7846918) q[0];
rz(-2.3214031) q[1];
sx q[1];
rz(-2.7892022) q[1];
sx q[1];
rz(-1.483451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950715) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(1.4875862) q[0];
rz(-pi) q[1];
rz(-2.5513068) q[2];
sx q[2];
rz(-1.4795408) q[2];
sx q[2];
rz(2.0412738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2918313) q[1];
sx q[1];
rz(-1.9007678) q[1];
sx q[1];
rz(-2.989819) q[1];
rz(2.5914707) q[3];
sx q[3];
rz(-1.2391483) q[3];
sx q[3];
rz(-1.5105997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(-1.8971987) q[2];
rz(2.7164298) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(-1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72220951) q[0];
sx q[0];
rz(-2.6477551) q[0];
sx q[0];
rz(-2.7710932) q[0];
rz(1.1100948) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.3409021) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6275218) q[0];
sx q[0];
rz(-2.7033269) q[0];
sx q[0];
rz(-2.1379495) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71435931) q[2];
sx q[2];
rz(-1.4673125) q[2];
sx q[2];
rz(-1.4504745) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95883137) q[1];
sx q[1];
rz(-1.8204096) q[1];
sx q[1];
rz(1.1968489) q[1];
rz(0.12124664) q[3];
sx q[3];
rz(-1.5910501) q[3];
sx q[3];
rz(0.95758394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5122539) q[2];
sx q[2];
rz(-2.9014355) q[2];
sx q[2];
rz(1.6188999) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.4168134) q[3];
sx q[3];
rz(-0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387977) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(2.8181656) q[1];
sx q[1];
rz(-1.2558116) q[1];
sx q[1];
rz(-1.5423923) q[1];
rz(1.1644438) q[2];
sx q[2];
rz(-1.5487557) q[2];
sx q[2];
rz(1.1974481) q[2];
rz(1.0988416) q[3];
sx q[3];
rz(-1.9190211) q[3];
sx q[3];
rz(2.2091051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
