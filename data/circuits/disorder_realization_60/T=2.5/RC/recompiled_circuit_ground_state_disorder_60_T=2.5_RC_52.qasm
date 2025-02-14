OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.056501) q[0];
sx q[0];
rz(-2.8488475) q[0];
sx q[0];
rz(-2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(-1.1408495) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5802683) q[0];
sx q[0];
rz(-1.394857) q[0];
sx q[0];
rz(-0.12138155) q[0];
rz(-pi) q[1];
rz(-1.8542641) q[2];
sx q[2];
rz(-1.8745443) q[2];
sx q[2];
rz(-1.2992573) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4881058) q[1];
sx q[1];
rz(-1.7035023) q[1];
sx q[1];
rz(-2.7567171) q[1];
x q[2];
rz(1.2841746) q[3];
sx q[3];
rz(-2.686073) q[3];
sx q[3];
rz(-2.974433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4484278) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(1.9654174) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(-2.6068408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.46928826) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(-0.09672673) q[1];
sx q[1];
rz(-1.6177142) q[1];
sx q[1];
rz(1.2299889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7589129) q[0];
sx q[0];
rz(-1.935674) q[0];
sx q[0];
rz(-1.7880102) q[0];
rz(-pi) q[1];
rz(-2.582091) q[2];
sx q[2];
rz(-1.5299262) q[2];
sx q[2];
rz(0.80383435) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7736771) q[1];
sx q[1];
rz(-1.8080447) q[1];
sx q[1];
rz(-1.3648454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1925063) q[3];
sx q[3];
rz(-1.6069901) q[3];
sx q[3];
rz(-0.32603273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35473287) q[2];
sx q[2];
rz(-0.41149461) q[2];
sx q[2];
rz(3.1301609) q[2];
rz(1.8276021) q[3];
sx q[3];
rz(-1.3649789) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74362022) q[0];
sx q[0];
rz(-0.13298661) q[0];
sx q[0];
rz(-0.61070329) q[0];
rz(0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(0.59649831) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4134408) q[0];
sx q[0];
rz(-0.96609173) q[0];
sx q[0];
rz(-0.20264969) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0192415) q[2];
sx q[2];
rz(-0.62659133) q[2];
sx q[2];
rz(-2.5950876) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8104659) q[1];
sx q[1];
rz(-2.2743723) q[1];
sx q[1];
rz(1.9845647) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5174808) q[3];
sx q[3];
rz(-2.5870596) q[3];
sx q[3];
rz(-2.7729101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67768031) q[2];
sx q[2];
rz(-1.3119421) q[2];
sx q[2];
rz(0.00027351969) q[2];
rz(1.4759981) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(-0.10703787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.45604712) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(2.1641425) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(-1.5553442) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.453131) q[0];
sx q[0];
rz(-1.5785839) q[0];
sx q[0];
rz(1.5611737) q[0];
rz(-0.78537099) q[2];
sx q[2];
rz(-1.6024688) q[2];
sx q[2];
rz(-0.23508628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11862893) q[1];
sx q[1];
rz(-1.5852155) q[1];
sx q[1];
rz(2.6447536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8951696) q[3];
sx q[3];
rz(-1.0895673) q[3];
sx q[3];
rz(-0.53477609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1996475) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(0.23435782) q[2];
rz(0.50218454) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(-1.3932047) q[0];
rz(0.69270095) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(2.0511973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4596953) q[0];
sx q[0];
rz(-2.5715264) q[0];
sx q[0];
rz(-1.1718114) q[0];
x q[1];
rz(-2.2545891) q[2];
sx q[2];
rz(-0.74216026) q[2];
sx q[2];
rz(-0.61486828) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0648374) q[1];
sx q[1];
rz(-1.7211815) q[1];
sx q[1];
rz(-1.9888179) q[1];
rz(1.0308068) q[3];
sx q[3];
rz(-1.705369) q[3];
sx q[3];
rz(-2.7887044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56110567) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(-2.6055276) q[2];
rz(0.29340336) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368161) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(0.099040898) q[0];
rz(-1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(-1.7031857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.563614) q[0];
sx q[0];
rz(-0.7398842) q[0];
sx q[0];
rz(-1.1881962) q[0];
rz(0.2834592) q[2];
sx q[2];
rz(-1.4999031) q[2];
sx q[2];
rz(-2.9170582) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21675303) q[1];
sx q[1];
rz(-2.3509563) q[1];
sx q[1];
rz(-0.032921493) q[1];
x q[2];
rz(-2.766633) q[3];
sx q[3];
rz(-0.31878933) q[3];
sx q[3];
rz(-2.6998883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1556039) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(-2.3206553) q[2];
rz(-1.692903) q[3];
sx q[3];
rz(-0.71805787) q[3];
sx q[3];
rz(-1.4887571) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(0.50265092) q[0];
rz(1.2333168) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(-0.9800235) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21979688) q[0];
sx q[0];
rz(-1.3338519) q[0];
sx q[0];
rz(-1.8264997) q[0];
rz(-pi) q[1];
rz(-0.77620164) q[2];
sx q[2];
rz(-1.7817678) q[2];
sx q[2];
rz(-2.7609776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7240643) q[1];
sx q[1];
rz(-1.8019946) q[1];
sx q[1];
rz(-0.75059143) q[1];
rz(1.3408889) q[3];
sx q[3];
rz(-2.3977444) q[3];
sx q[3];
rz(0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9194453) q[2];
sx q[2];
rz(-1.28747) q[2];
sx q[2];
rz(1.3745098) q[2];
rz(-1.8253271) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.114349) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(0.91841665) q[0];
rz(0.20485993) q[1];
sx q[1];
rz(-2.0826976) q[1];
sx q[1];
rz(-1.6580261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99824725) q[0];
sx q[0];
rz(-1.7172617) q[0];
sx q[0];
rz(1.2781758) q[0];
rz(-pi) q[1];
rz(-0.81088938) q[2];
sx q[2];
rz(-1.4810908) q[2];
sx q[2];
rz(2.9765338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10459722) q[1];
sx q[1];
rz(-0.87338398) q[1];
sx q[1];
rz(-0.11499494) q[1];
rz(-1.5213883) q[3];
sx q[3];
rz(-0.8522343) q[3];
sx q[3];
rz(-0.75203943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1398937) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(-1.3341382) q[2];
rz(-2.259518) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(-1.849256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.081721574) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(-3.1234142) q[0];
rz(2.7320618) q[1];
sx q[1];
rz(-1.4298226) q[1];
sx q[1];
rz(3.0407564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5715953) q[0];
sx q[0];
rz(-2.2048108) q[0];
sx q[0];
rz(1.5279395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11282044) q[2];
sx q[2];
rz(-2.718528) q[2];
sx q[2];
rz(-0.12020883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68647766) q[1];
sx q[1];
rz(-1.9338738) q[1];
sx q[1];
rz(-2.8836714) q[1];
rz(-pi) q[2];
rz(-1.8180679) q[3];
sx q[3];
rz(-1.4535289) q[3];
sx q[3];
rz(-0.7922678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
sx q[2];
rz(1.8936554) q[2];
rz(-1.2248056) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(-3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.4757272) q[0];
sx q[0];
rz(-2.8210848) q[0];
rz(0.39101741) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(1.5919707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3452197) q[0];
sx q[0];
rz(-1.380305) q[0];
sx q[0];
rz(2.9389771) q[0];
rz(-pi) q[1];
rz(-2.8587927) q[2];
sx q[2];
rz(-0.53164266) q[2];
sx q[2];
rz(-0.34257364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84234778) q[1];
sx q[1];
rz(-1.4834373) q[1];
sx q[1];
rz(-1.9398938) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1248766) q[3];
sx q[3];
rz(-2.7191945) q[3];
sx q[3];
rz(2.0579541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0512507) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-0.40038294) q[2];
rz(-0.23614899) q[3];
sx q[3];
rz(-3.0381687) q[3];
sx q[3];
rz(2.1307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-2.6564468) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(-1.0322124) q[2];
sx q[2];
rz(-2.9001424) q[2];
sx q[2];
rz(2.7912414) q[2];
rz(-1.7988206) q[3];
sx q[3];
rz(-1.6327471) q[3];
sx q[3];
rz(-3.0625797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
