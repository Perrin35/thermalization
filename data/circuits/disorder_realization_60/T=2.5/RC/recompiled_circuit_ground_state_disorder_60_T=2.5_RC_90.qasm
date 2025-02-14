OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(2.2226287) q[0];
rz(-0.38209823) q[1];
sx q[1];
rz(-2.9736019) q[1];
sx q[1];
rz(1.1408495) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1699693) q[0];
sx q[0];
rz(-0.21339082) q[0];
sx q[0];
rz(2.1687228) q[0];
x q[1];
rz(2.4130526) q[2];
sx q[2];
rz(-2.72914) q[2];
sx q[2];
rz(-1.0701239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.39841043) q[1];
sx q[1];
rz(-2.7355621) q[1];
sx q[1];
rz(2.7999785) q[1];
rz(-1.857418) q[3];
sx q[3];
rz(-2.686073) q[3];
sx q[3];
rz(0.16715967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69316489) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(-1.9654174) q[2];
rz(3.0991992) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(-0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46928826) q[0];
sx q[0];
rz(-2.7870218) q[0];
sx q[0];
rz(2.9571423) q[0];
rz(3.0448659) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(1.9116037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.205015) q[0];
sx q[0];
rz(-0.42213556) q[0];
sx q[0];
rz(-0.5136712) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0646992) q[2];
sx q[2];
rz(-0.56083365) q[2];
sx q[2];
rz(-0.70181134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0470386) q[1];
sx q[1];
rz(-2.8287005) q[1];
sx q[1];
rz(-0.70200664) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6685247) q[3];
sx q[3];
rz(-2.7616581) q[3];
sx q[3];
rz(-1.987628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979724) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(-0.61070329) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(2.5450943) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72815182) q[0];
sx q[0];
rz(-2.1755009) q[0];
sx q[0];
rz(-0.20264969) q[0];
x q[1];
rz(-2.1488068) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(1.3959194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92745078) q[1];
sx q[1];
rz(-2.3436556) q[1];
sx q[1];
rz(-0.44254704) q[1];
x q[2];
rz(1.6241118) q[3];
sx q[3];
rz(-0.55453306) q[3];
sx q[3];
rz(-0.36868251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(3.1413191) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(0.10703787) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(-2.1641425) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(1.5862484) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.259183) q[0];
sx q[0];
rz(-1.5611739) q[0];
sx q[0];
rz(3.1338047) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6155717) q[2];
sx q[2];
rz(-2.3556659) q[2];
sx q[2];
rz(1.3040552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4787538) q[1];
sx q[1];
rz(-0.49703076) q[1];
sx q[1];
rz(-0.030244002) q[1];
rz(-pi) q[2];
rz(2.0078251) q[3];
sx q[3];
rz(-0.53621447) q[3];
sx q[3];
rz(1.0325583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.9419452) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-2.9072348) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-1.0415223) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(1.748388) q[0];
rz(-0.69270095) q[1];
sx q[1];
rz(-2.0557978) q[1];
sx q[1];
rz(1.0903953) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22995725) q[0];
sx q[0];
rz(-1.3595694) q[0];
sx q[0];
rz(-1.0372355) q[0];
rz(-pi) q[1];
rz(0.95280052) q[2];
sx q[2];
rz(-2.0119466) q[2];
sx q[2];
rz(-1.49681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0648374) q[1];
sx q[1];
rz(-1.4204111) q[1];
sx q[1];
rz(1.1527747) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15655915) q[3];
sx q[3];
rz(-1.0362175) q[3];
sx q[3];
rz(-1.8434356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(2.6055276) q[2];
rz(2.8481893) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(-2.6611879) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90477657) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(3.0425518) q[0];
rz(-2.1095443) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(-1.438407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0625298) q[0];
sx q[0];
rz(-2.2465116) q[0];
sx q[0];
rz(-2.8131301) q[0];
x q[1];
rz(2.8929404) q[2];
sx q[2];
rz(-2.8496345) q[2];
sx q[2];
rz(-1.1077293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.878053) q[1];
sx q[1];
rz(-2.3608853) q[1];
sx q[1];
rz(1.6040463) q[1];
rz(2.766633) q[3];
sx q[3];
rz(-0.31878933) q[3];
sx q[3];
rz(2.6998883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1556039) q[2];
sx q[2];
rz(-0.5039379) q[2];
sx q[2];
rz(-0.8209374) q[2];
rz(1.692903) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(0.50265092) q[0];
rz(-1.9082759) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(2.1615692) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7293046) q[0];
sx q[0];
rz(-1.8192023) q[0];
sx q[0];
rz(-0.24459837) q[0];
x q[1];
rz(-2.8449242) q[2];
sx q[2];
rz(-0.7985332) q[2];
sx q[2];
rz(-0.980033) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.087638559) q[1];
sx q[1];
rz(-0.77869895) q[1];
sx q[1];
rz(-0.33233541) q[1];
rz(1.8007038) q[3];
sx q[3];
rz(-0.74384825) q[3];
sx q[3];
rz(0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.7670828) q[2];
rz(-1.8253271) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(-2.159806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027243622) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(-0.91841665) q[0];
rz(-2.9367327) q[1];
sx q[1];
rz(-2.0826976) q[1];
sx q[1];
rz(1.4835666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0201538) q[0];
sx q[0];
rz(-2.8153072) q[0];
sx q[0];
rz(2.043528) q[0];
rz(1.4409562) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(1.4997327) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0682448) q[1];
sx q[1];
rz(-2.4363344) q[1];
sx q[1];
rz(1.4346992) q[1];
rz(-0.056413944) q[3];
sx q[3];
rz(-2.4216363) q[3];
sx q[3];
rz(0.82701339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1398937) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(-1.8074544) q[2];
rz(2.259518) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(1.2923366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598711) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(0.01817848) q[0];
rz(0.40953088) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(-0.10083625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4977328) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(-3.0833901) q[0];
x q[1];
rz(1.5201496) q[2];
sx q[2];
rz(-1.9910004) q[2];
sx q[2];
rz(0.24383185) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3253357) q[1];
sx q[1];
rz(-0.44202572) q[1];
sx q[1];
rz(-2.1620552) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0193897) q[3];
sx q[3];
rz(-2.8684385) q[3];
sx q[3];
rz(2.7970683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8044901) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(-1.2479372) q[2];
rz(1.2248056) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(2.8210848) q[0];
rz(2.7505752) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(-1.5919707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79637291) q[0];
sx q[0];
rz(-1.7612877) q[0];
sx q[0];
rz(0.20261554) q[0];
rz(-pi) q[1];
rz(1.7334601) q[2];
sx q[2];
rz(-1.062359) q[2];
sx q[2];
rz(0.66772738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.446881) q[1];
sx q[1];
rz(-1.9384192) q[1];
sx q[1];
rz(3.047961) q[1];
x q[2];
rz(0.42234588) q[3];
sx q[3];
rz(-1.5776488) q[3];
sx q[3];
rz(0.50240483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.090342) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-2.7412097) q[2];
rz(-0.23614899) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(1.0108488) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26265963) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-2.6564468) q[1];
sx q[1];
rz(-2.3427675) q[1];
sx q[1];
rz(2.6642703) q[1];
rz(1.7791228) q[2];
sx q[2];
rz(-1.6937509) q[2];
sx q[2];
rz(-1.3954336) q[2];
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
