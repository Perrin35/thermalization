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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(2.5757134) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(-1.3574418) q[1];
sx q[1];
rz(0.53425962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256828) q[0];
sx q[0];
rz(-0.74640154) q[0];
sx q[0];
rz(-1.9749667) q[0];
rz(-pi) q[1];
rz(1.8741605) q[2];
sx q[2];
rz(-2.144226) q[2];
sx q[2];
rz(-1.0614741) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6837343) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(0.079816743) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47655063) q[3];
sx q[3];
rz(-2.635987) q[3];
sx q[3];
rz(-1.4534392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.013022097) q[2];
sx q[2];
rz(-0.2505005) q[2];
sx q[2];
rz(-2.8924083) q[2];
rz(1.7976044) q[3];
sx q[3];
rz(-1.5430217) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473815) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(2.1517854) q[0];
rz(-0.79065943) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(2.8299423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79651901) q[0];
sx q[0];
rz(-0.88788827) q[0];
sx q[0];
rz(0.89181283) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0333854) q[2];
sx q[2];
rz(-1.3851067) q[2];
sx q[2];
rz(-0.78150392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5492925) q[1];
sx q[1];
rz(-2.8452497) q[1];
sx q[1];
rz(-2.1974988) q[1];
x q[2];
rz(0.95590638) q[3];
sx q[3];
rz(-1.8528189) q[3];
sx q[3];
rz(-0.36263719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(0.95138335) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15762873) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(0.11446318) q[0];
rz(1.0022256) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(-0.1618298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1755668) q[0];
sx q[0];
rz(-1.1822261) q[0];
sx q[0];
rz(2.2184664) q[0];
x q[1];
rz(1.6587018) q[2];
sx q[2];
rz(-1.958333) q[2];
sx q[2];
rz(-1.0532925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.964965) q[1];
sx q[1];
rz(-1.8102268) q[1];
sx q[1];
rz(2.5188732) q[1];
rz(-pi) q[2];
rz(0.39997081) q[3];
sx q[3];
rz(-2.2999527) q[3];
sx q[3];
rz(0.10310452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3942922) q[2];
sx q[2];
rz(-0.44555274) q[2];
sx q[2];
rz(-3.0008345) q[2];
rz(0.55177871) q[3];
sx q[3];
rz(-1.2740302) q[3];
sx q[3];
rz(0.75132918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20525876) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(-0.67489135) q[0];
rz(-0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-0.45796576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5602011) q[0];
sx q[0];
rz(-1.0835674) q[0];
sx q[0];
rz(-2.013252) q[0];
rz(-pi) q[1];
rz(-2.8343924) q[2];
sx q[2];
rz(-1.1591025) q[2];
sx q[2];
rz(1.9780985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1110036) q[1];
sx q[1];
rz(-1.9408424) q[1];
sx q[1];
rz(1.7683328) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9513957) q[3];
sx q[3];
rz(-0.51589291) q[3];
sx q[3];
rz(-2.7864393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13350479) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(-1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.1595681) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(1.9200448) q[0];
rz(0.23280652) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(2.1585042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670416) q[0];
sx q[0];
rz(-2.3877151) q[0];
sx q[0];
rz(2.3874832) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7295425) q[2];
sx q[2];
rz(-0.6266784) q[2];
sx q[2];
rz(1.145663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5447588) q[1];
sx q[1];
rz(-1.7167517) q[1];
sx q[1];
rz(-2.9533141) q[1];
rz(-1.9931204) q[3];
sx q[3];
rz(-1.4185126) q[3];
sx q[3];
rz(1.9773169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7588707) q[2];
sx q[2];
rz(-1.4843586) q[2];
sx q[2];
rz(-2.7663084) q[2];
rz(-1.075607) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4009092) q[0];
sx q[0];
rz(-1.704819) q[0];
sx q[0];
rz(1.4373454) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(-2.5934503) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63614908) q[0];
sx q[0];
rz(-0.6983122) q[0];
sx q[0];
rz(2.8092119) q[0];
rz(-pi) q[1];
rz(-2.5896543) q[2];
sx q[2];
rz(-1.9935529) q[2];
sx q[2];
rz(1.1814337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9029451) q[1];
sx q[1];
rz(-2.3460707) q[1];
sx q[1];
rz(-2.5832157) q[1];
x q[2];
rz(-2.4924566) q[3];
sx q[3];
rz(-0.86926354) q[3];
sx q[3];
rz(2.8809406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6846201) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(-0.6753298) q[2];
rz(-1.5843377) q[3];
sx q[3];
rz(-1.3074338) q[3];
sx q[3];
rz(-0.61711446) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5439344) q[0];
sx q[0];
rz(-1.958853) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(2.0248263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1710715) q[0];
sx q[0];
rz(-1.9200033) q[0];
sx q[0];
rz(1.1961351) q[0];
x q[1];
rz(0.26979763) q[2];
sx q[2];
rz(-1.0493663) q[2];
sx q[2];
rz(0.21873378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6903444) q[1];
sx q[1];
rz(-0.4745698) q[1];
sx q[1];
rz(0.75466538) q[1];
rz(-pi) q[2];
rz(-2.6855822) q[3];
sx q[3];
rz(-1.6263762) q[3];
sx q[3];
rz(-1.0694651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(0.35487077) q[2];
rz(-0.76534671) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(-0.089546831) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(0.76559693) q[0];
rz(-0.87144026) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(1.3227468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3324748) q[0];
sx q[0];
rz(-1.6226543) q[0];
sx q[0];
rz(2.8807667) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23643146) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(-2.1588391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31963667) q[1];
sx q[1];
rz(-1.0624891) q[1];
sx q[1];
rz(-1.5519358) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20669528) q[3];
sx q[3];
rz(-1.350025) q[3];
sx q[3];
rz(-1.4109991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(0.22171177) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(-2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0677277) q[0];
sx q[0];
rz(-1.2465957) q[0];
sx q[0];
rz(2.8644323) q[0];
rz(-0.74835888) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(-1.1433196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4244014) q[0];
sx q[0];
rz(-0.36866185) q[0];
sx q[0];
rz(2.1564846) q[0];
x q[1];
rz(-0.21094811) q[2];
sx q[2];
rz(-2.2648513) q[2];
sx q[2];
rz(-1.2488493) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6618234) q[1];
sx q[1];
rz(-2.7854438) q[1];
sx q[1];
rz(-1.1060064) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42871126) q[3];
sx q[3];
rz(-0.39434163) q[3];
sx q[3];
rz(-1.1129781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9144168) q[2];
sx q[2];
rz(-2.1180426) q[2];
sx q[2];
rz(0.048967036) q[2];
rz(-1.0660727) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(-0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-0.019512026) q[0];
rz(-2.3628269) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(-1.5628409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21694788) q[0];
sx q[0];
rz(-2.9836982) q[0];
sx q[0];
rz(2.1720706) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1807272) q[2];
sx q[2];
rz(-1.4205407) q[2];
sx q[2];
rz(-0.021266887) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9924888) q[1];
sx q[1];
rz(-1.1380151) q[1];
sx q[1];
rz(-0.91723587) q[1];
x q[2];
rz(2.9928777) q[3];
sx q[3];
rz(-0.84245719) q[3];
sx q[3];
rz(3.0535798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(0.16555244) q[2];
rz(-0.60756573) q[3];
sx q[3];
rz(-1.8568042) q[3];
sx q[3];
rz(0.46835381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-1.7623357) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
rz(1.5070076) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(2.9956837) q[2];
sx q[2];
rz(-0.78073954) q[2];
sx q[2];
rz(-2.2013046) q[2];
rz(-0.60579469) q[3];
sx q[3];
rz(-1.4069362) q[3];
sx q[3];
rz(-2.8154472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
