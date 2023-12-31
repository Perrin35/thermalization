OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(-1.5024827) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26412548) q[0];
sx q[0];
rz(-1.2116417) q[0];
sx q[0];
rz(-2.9496664) q[0];
rz(-1.5760954) q[2];
sx q[2];
rz(-2.1359348) q[2];
sx q[2];
rz(2.0084755) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2734387) q[1];
sx q[1];
rz(-1.3248797) q[1];
sx q[1];
rz(-1.6383365) q[1];
rz(-pi) q[2];
rz(-1.1611657) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(-0.15256552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(0.65594977) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(-0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(2.7080652) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1996961) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.9641563) q[0];
rz(-pi) q[1];
rz(-0.29351182) q[2];
sx q[2];
rz(-1.8732757) q[2];
sx q[2];
rz(0.40700618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.05939535) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(-2.5897964) q[1];
x q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.1840608) q[3];
sx q[3];
rz(2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5269512) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(0.69765222) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.7720222) q[0];
rz(-1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(2.8799768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152552) q[0];
sx q[0];
rz(-1.5528423) q[0];
sx q[0];
rz(0.0096339027) q[0];
rz(-3.0171266) q[2];
sx q[2];
rz(-2.3111768) q[2];
sx q[2];
rz(0.24701961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9237823) q[1];
sx q[1];
rz(-0.79499704) q[1];
sx q[1];
rz(-1.4008821) q[1];
x q[2];
rz(1.5393799) q[3];
sx q[3];
rz(-2.0835702) q[3];
sx q[3];
rz(-0.98785066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(-0.20351163) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(-1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.2483695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49444775) q[0];
sx q[0];
rz(-1.6361423) q[0];
sx q[0];
rz(1.4739743) q[0];
x q[1];
rz(-2.4291971) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(1.4272387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(-0.65472366) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4183211) q[3];
sx q[3];
rz(-1.7687106) q[3];
sx q[3];
rz(-2.7151782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(2.3045585) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-2.2391589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9426614) q[0];
sx q[0];
rz(-2.5006602) q[0];
sx q[0];
rz(0.076365691) q[0];
rz(-pi) q[1];
rz(-2.6685153) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(-0.15854533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8946998) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(2.488399) q[1];
x q[2];
rz(-1.8875214) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9465785) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.1966594) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3301795) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(2.6384171) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(0.20176372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70226442) q[0];
sx q[0];
rz(-1.6730671) q[0];
sx q[0];
rz(-0.066145397) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1226419) q[2];
sx q[2];
rz(-1.1345703) q[2];
sx q[2];
rz(3.089038) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4469874) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(-2.6161731) q[1];
rz(-pi) q[2];
rz(2.0353349) q[3];
sx q[3];
rz(-0.89336508) q[3];
sx q[3];
rz(2.5708452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(2.8743437) q[2];
rz(0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(-2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(3.0840432) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(2.1988791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7563815) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-2.5368607) q[0];
x q[1];
rz(-0.16072388) q[2];
sx q[2];
rz(-1.7956927) q[2];
sx q[2];
rz(-2.7149534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49589866) q[1];
sx q[1];
rz(-1.2518479) q[1];
sx q[1];
rz(1.1369399) q[1];
rz(2.8816678) q[3];
sx q[3];
rz(-2.0689031) q[3];
sx q[3];
rz(-2.2181524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(-2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(-0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(2.8576635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89284183) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(0.52389223) q[0];
rz(-1.7173041) q[2];
sx q[2];
rz(-0.97587817) q[2];
sx q[2];
rz(-0.3301687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24208454) q[1];
sx q[1];
rz(-0.71166066) q[1];
sx q[1];
rz(-2.3942024) q[1];
rz(-2.4817804) q[3];
sx q[3];
rz(-1.5944578) q[3];
sx q[3];
rz(-0.82994474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.4607666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979249) q[0];
sx q[0];
rz(-2.857508) q[0];
sx q[0];
rz(-0.062505917) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.1968687) q[2];
sx q[2];
rz(1.1073081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31635346) q[1];
sx q[1];
rz(-1.5484527) q[1];
sx q[1];
rz(0.52492001) q[1];
rz(-pi) q[2];
rz(0.78482307) q[3];
sx q[3];
rz(-1.8514957) q[3];
sx q[3];
rz(-0.37898889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.573695) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(0.38063231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779553) q[0];
sx q[0];
rz(-1.5008238) q[0];
sx q[0];
rz(-1.3689343) q[0];
x q[1];
rz(0.62217525) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(3.0058793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3818647) q[1];
sx q[1];
rz(-1.7654395) q[1];
sx q[1];
rz(-1.0641644) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5467039) q[3];
sx q[3];
rz(-0.66685646) q[3];
sx q[3];
rz(0.24645933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0011065817) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(-0.040955695) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363591) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(1.8006051) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(-0.20053486) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
