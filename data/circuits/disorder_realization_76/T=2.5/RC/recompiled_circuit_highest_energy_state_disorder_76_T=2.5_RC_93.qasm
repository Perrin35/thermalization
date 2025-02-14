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
rz(-2.2657356) q[0];
sx q[0];
rz(-2.1729204) q[0];
sx q[0];
rz(0.92585603) q[0];
rz(-2.5246188) q[1];
sx q[1];
rz(-2.4763835) q[1];
sx q[1];
rz(1.8170504) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449555) q[0];
sx q[0];
rz(-0.60316604) q[0];
sx q[0];
rz(-0.11267333) q[0];
rz(2.994903) q[2];
sx q[2];
rz(-2.5518199) q[2];
sx q[2];
rz(-1.9751994) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.032925978) q[1];
sx q[1];
rz(-0.98331988) q[1];
sx q[1];
rz(1.0125748) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67456986) q[3];
sx q[3];
rz(-2.2534182) q[3];
sx q[3];
rz(-2.5404524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8806184) q[2];
sx q[2];
rz(-1.7558492) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51139128) q[0];
sx q[0];
rz(-1.8237317) q[0];
sx q[0];
rz(2.200101) q[0];
rz(-0.9990274) q[1];
sx q[1];
rz(-2.2143054) q[1];
sx q[1];
rz(-0.082854465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1076528) q[0];
sx q[0];
rz(-1.4705188) q[0];
sx q[0];
rz(2.4058008) q[0];
rz(-pi) q[1];
rz(-2.8993334) q[2];
sx q[2];
rz(-1.2238811) q[2];
sx q[2];
rz(1.0936979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1260316) q[1];
sx q[1];
rz(-1.120561) q[1];
sx q[1];
rz(-0.57029389) q[1];
rz(-2.4965246) q[3];
sx q[3];
rz(-1.0655902) q[3];
sx q[3];
rz(-0.92407214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2443709) q[2];
sx q[2];
rz(-1.3542078) q[2];
sx q[2];
rz(1.2574035) q[2];
rz(-2.2952378) q[3];
sx q[3];
rz(-2.9823163) q[3];
sx q[3];
rz(-2.4634821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427247) q[0];
sx q[0];
rz(-1.4582448) q[0];
sx q[0];
rz(0.22920907) q[0];
rz(0.21408679) q[1];
sx q[1];
rz(-2.2702718) q[1];
sx q[1];
rz(0.3814989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70951916) q[0];
sx q[0];
rz(-0.95147485) q[0];
sx q[0];
rz(1.1362057) q[0];
rz(1.7842641) q[2];
sx q[2];
rz(-0.74356438) q[2];
sx q[2];
rz(1.5876169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.72813624) q[1];
sx q[1];
rz(-1.6237139) q[1];
sx q[1];
rz(0.028832988) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1274453) q[3];
sx q[3];
rz(-0.34711829) q[3];
sx q[3];
rz(1.0747758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7233589) q[2];
sx q[2];
rz(-0.25622076) q[2];
sx q[2];
rz(2.5733433) q[2];
rz(1.7226284) q[3];
sx q[3];
rz(-1.7122372) q[3];
sx q[3];
rz(1.0770146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028246183) q[0];
sx q[0];
rz(-1.1320817) q[0];
sx q[0];
rz(-2.8908253) q[0];
rz(-3.057042) q[1];
sx q[1];
rz(-2.0471579) q[1];
sx q[1];
rz(1.9409404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687209) q[0];
sx q[0];
rz(-0.78470147) q[0];
sx q[0];
rz(1.0502104) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0122288) q[2];
sx q[2];
rz(-1.7953292) q[2];
sx q[2];
rz(0.44432902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29087092) q[1];
sx q[1];
rz(-2.4230036) q[1];
sx q[1];
rz(0.1512089) q[1];
rz(2.0968998) q[3];
sx q[3];
rz(-1.8013489) q[3];
sx q[3];
rz(1.5073204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4667929) q[2];
sx q[2];
rz(-2.0189221) q[2];
sx q[2];
rz(1.8505081) q[2];
rz(-2.71991) q[3];
sx q[3];
rz(-0.45160523) q[3];
sx q[3];
rz(1.565056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6186433) q[0];
sx q[0];
rz(-0.72294253) q[0];
sx q[0];
rz(-1.6408386) q[0];
rz(0.067642637) q[1];
sx q[1];
rz(-1.5875971) q[1];
sx q[1];
rz(1.1281475) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2826506) q[0];
sx q[0];
rz(-2.9500486) q[0];
sx q[0];
rz(-2.9454977) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7893547) q[2];
sx q[2];
rz(-2.0697429) q[2];
sx q[2];
rz(-1.3832472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6349244) q[1];
sx q[1];
rz(-1.4574058) q[1];
sx q[1];
rz(-0.054111295) q[1];
x q[2];
rz(1.2913843) q[3];
sx q[3];
rz(-1.3456723) q[3];
sx q[3];
rz(-2.1373419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1077914) q[2];
sx q[2];
rz(-2.0267603) q[2];
sx q[2];
rz(2.4647253) q[2];
rz(-0.90773165) q[3];
sx q[3];
rz(-1.9151442) q[3];
sx q[3];
rz(-0.41456732) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.2409302) q[1];
sx q[1];
rz(-2.2878094) q[1];
sx q[1];
rz(0.21533899) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8148635) q[0];
sx q[0];
rz(-1.0007326) q[0];
sx q[0];
rz(0.22320052) q[0];
rz(-pi) q[1];
rz(-1.239993) q[2];
sx q[2];
rz(-1.5546521) q[2];
sx q[2];
rz(0.46260297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78241888) q[1];
sx q[1];
rz(-2.4451588) q[1];
sx q[1];
rz(2.7553808) q[1];
rz(-1.9086544) q[3];
sx q[3];
rz(-0.90150276) q[3];
sx q[3];
rz(-0.026642628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.874959) q[2];
sx q[2];
rz(-1.6338438) q[2];
sx q[2];
rz(-0.61666644) q[2];
rz(-2.941361) q[3];
sx q[3];
rz(-2.4112406) q[3];
sx q[3];
rz(2.0088137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(3.1277593) q[0];
sx q[0];
rz(-0.56202373) q[0];
sx q[0];
rz(-0.01509893) q[0];
rz(3.0191782) q[1];
sx q[1];
rz(-1.2950803) q[1];
sx q[1];
rz(1.0282358) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434179) q[0];
sx q[0];
rz(-2.2255996) q[0];
sx q[0];
rz(-1.6408987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2790658) q[2];
sx q[2];
rz(-2.5030067) q[2];
sx q[2];
rz(-1.2995468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15785698) q[1];
sx q[1];
rz(-0.4326371) q[1];
sx q[1];
rz(-0.067072596) q[1];
rz(0.7042775) q[3];
sx q[3];
rz(-0.40474328) q[3];
sx q[3];
rz(1.2621439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.60468173) q[2];
sx q[2];
rz(-1.2102419) q[2];
sx q[2];
rz(-0.49986419) q[2];
rz(2.8258421) q[3];
sx q[3];
rz(-0.67086589) q[3];
sx q[3];
rz(1.0189112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4278118) q[0];
sx q[0];
rz(-1.2232057) q[0];
sx q[0];
rz(-2.1977303) q[0];
rz(-1.7367412) q[1];
sx q[1];
rz(-1.8753139) q[1];
sx q[1];
rz(-2.2199383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93225828) q[0];
sx q[0];
rz(-2.8186322) q[0];
sx q[0];
rz(0.17318101) q[0];
x q[1];
rz(-1.2723075) q[2];
sx q[2];
rz(-1.8137852) q[2];
sx q[2];
rz(1.1557494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6347305) q[1];
sx q[1];
rz(-0.94351879) q[1];
sx q[1];
rz(-0.92232134) q[1];
x q[2];
rz(-1.0579487) q[3];
sx q[3];
rz(-1.5745224) q[3];
sx q[3];
rz(1.8029574) q[3];
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
rz(-1.5550295) q[2];
rz(-1.0541213) q[3];
sx q[3];
rz(-1.328822) q[3];
sx q[3];
rz(-1.4012977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15040511) q[0];
sx q[0];
rz(-2.0441971) q[0];
sx q[0];
rz(0.64252585) q[0];
rz(1.7036899) q[1];
sx q[1];
rz(-0.75690401) q[1];
sx q[1];
rz(-2.9978851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3966345) q[0];
sx q[0];
rz(-1.8734583) q[0];
sx q[0];
rz(-2.55654) q[0];
x q[1];
rz(2.2177198) q[2];
sx q[2];
rz(-2.0028186) q[2];
sx q[2];
rz(0.96521711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2255114) q[1];
sx q[1];
rz(-0.71174946) q[1];
sx q[1];
rz(2.7307778) q[1];
rz(0.71227422) q[3];
sx q[3];
rz(-0.55432075) q[3];
sx q[3];
rz(-1.2620776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5370499) q[2];
sx q[2];
rz(-1.9436676) q[2];
sx q[2];
rz(2.878888) q[2];
rz(2.883834) q[3];
sx q[3];
rz(-2.7596605) q[3];
sx q[3];
rz(2.2885382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2951374) q[0];
sx q[0];
rz(-1.6520123) q[0];
sx q[0];
rz(1.9211796) q[0];
rz(-2.659761) q[1];
sx q[1];
rz(-2.316663) q[1];
sx q[1];
rz(0.89734546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2663956) q[0];
sx q[0];
rz(-1.0561854) q[0];
sx q[0];
rz(-1.4304016) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32081066) q[2];
sx q[2];
rz(-1.6751373) q[2];
sx q[2];
rz(3.0719245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1188965) q[1];
sx q[1];
rz(-1.6204981) q[1];
sx q[1];
rz(-1.3124052) q[1];
rz(2.1796024) q[3];
sx q[3];
rz(-1.4856087) q[3];
sx q[3];
rz(0.60768581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64409488) q[2];
sx q[2];
rz(-0.92659014) q[2];
sx q[2];
rz(-3.0268055) q[2];
rz(-2.3980906) q[3];
sx q[3];
rz(-1.2657575) q[3];
sx q[3];
rz(0.44873294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262309) q[0];
sx q[0];
rz(-2.3804433) q[0];
sx q[0];
rz(-0.95809715) q[0];
rz(-1.6356946) q[1];
sx q[1];
rz(-2.0151357) q[1];
sx q[1];
rz(-1.9912079) q[1];
rz(2.4534879) q[2];
sx q[2];
rz(-2.2672014) q[2];
sx q[2];
rz(-1.3939569) q[2];
rz(-1.7382449) q[3];
sx q[3];
rz(-1.9994152) q[3];
sx q[3];
rz(-2.216702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
