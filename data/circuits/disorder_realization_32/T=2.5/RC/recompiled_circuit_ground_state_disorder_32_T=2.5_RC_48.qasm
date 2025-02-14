OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(7.5543348) q[0];
sx q[0];
rz(13.332097) q[0];
rz(2.7453121) q[1];
sx q[1];
rz(-3.0488465) q[1];
sx q[1];
rz(2.588811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6792949) q[0];
sx q[0];
rz(-0.57749417) q[0];
sx q[0];
rz(2.3403843) q[0];
rz(0.80443126) q[2];
sx q[2];
rz(-1.0519769) q[2];
sx q[2];
rz(0.08229736) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0243135) q[1];
sx q[1];
rz(-1.8872042) q[1];
sx q[1];
rz(-0.98543075) q[1];
rz(0.90867282) q[3];
sx q[3];
rz(-1.931802) q[3];
sx q[3];
rz(-0.26546016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3679686) q[2];
sx q[2];
rz(-1.6398733) q[2];
sx q[2];
rz(0.088851301) q[2];
rz(0.39696524) q[3];
sx q[3];
rz(-2.1019955) q[3];
sx q[3];
rz(-1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6615768) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(2.2870824) q[0];
rz(2.5201733) q[1];
sx q[1];
rz(-0.3365376) q[1];
sx q[1];
rz(1.052676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7490279) q[0];
sx q[0];
rz(-1.6727007) q[0];
sx q[0];
rz(0.13691957) q[0];
rz(-pi) q[1];
rz(2.7268953) q[2];
sx q[2];
rz(-1.9837579) q[2];
sx q[2];
rz(-0.1916445) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32089511) q[1];
sx q[1];
rz(-1.7716494) q[1];
sx q[1];
rz(-0.83495514) q[1];
x q[2];
rz(-0.0016453513) q[3];
sx q[3];
rz(-2.1747394) q[3];
sx q[3];
rz(1.2206447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66775995) q[2];
sx q[2];
rz(-2.2809873) q[2];
sx q[2];
rz(-2.1916126) q[2];
rz(-2.1330323) q[3];
sx q[3];
rz(-0.63298321) q[3];
sx q[3];
rz(2.8620201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94473332) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(-1.3470294) q[0];
rz(2.0836209) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(-0.11014858) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5025217) q[0];
sx q[0];
rz(-0.57146954) q[0];
sx q[0];
rz(-2.1854464) q[0];
x q[1];
rz(3.0579849) q[2];
sx q[2];
rz(-1.7183398) q[2];
sx q[2];
rz(1.7718441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5837142) q[1];
sx q[1];
rz(-1.4643351) q[1];
sx q[1];
rz(0.6006247) q[1];
rz(-pi) q[2];
rz(-0.83602704) q[3];
sx q[3];
rz(-0.6202094) q[3];
sx q[3];
rz(2.7171752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9614253) q[2];
sx q[2];
rz(-1.5632997) q[2];
sx q[2];
rz(0.041672826) q[2];
rz(-2.9336119) q[3];
sx q[3];
rz(-0.22170034) q[3];
sx q[3];
rz(2.8878133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5792849) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(-2.1606309) q[0];
rz(0.43308577) q[1];
sx q[1];
rz(-1.667645) q[1];
sx q[1];
rz(-1.513419) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53674066) q[0];
sx q[0];
rz(-1.3910595) q[0];
sx q[0];
rz(0.40412235) q[0];
rz(-pi) q[1];
rz(3.1279706) q[2];
sx q[2];
rz(-0.8994461) q[2];
sx q[2];
rz(-1.4289795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1867552) q[1];
sx q[1];
rz(-1.4680956) q[1];
sx q[1];
rz(2.0689194) q[1];
rz(-pi) q[2];
rz(1.9497197) q[3];
sx q[3];
rz(-0.66693587) q[3];
sx q[3];
rz(2.2569424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.607434) q[2];
sx q[2];
rz(-2.6223493) q[2];
sx q[2];
rz(0.81238166) q[2];
rz(1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(-1.2468503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3355376) q[0];
sx q[0];
rz(-1.0220818) q[0];
sx q[0];
rz(1.442765) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.6280326) q[1];
sx q[1];
rz(-2.3853669) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4837912) q[0];
sx q[0];
rz(-1.4157802) q[0];
sx q[0];
rz(-2.20138) q[0];
rz(-pi) q[1];
rz(2.4632023) q[2];
sx q[2];
rz(-2.3783461) q[2];
sx q[2];
rz(0.89707546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5797983) q[1];
sx q[1];
rz(-1.4358014) q[1];
sx q[1];
rz(-1.4827824) q[1];
x q[2];
rz(1.9147768) q[3];
sx q[3];
rz(-1.5369253) q[3];
sx q[3];
rz(2.276941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68236399) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(0.96561042) q[2];
rz(-2.5718578) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(1.6857356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4023034) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(1.385561) q[0];
rz(2.3784474) q[1];
sx q[1];
rz(-1.6940073) q[1];
sx q[1];
rz(-0.10173434) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4674326) q[0];
sx q[0];
rz(-1.6629574) q[0];
sx q[0];
rz(2.1919495) q[0];
rz(2.7587682) q[2];
sx q[2];
rz(-1.5271795) q[2];
sx q[2];
rz(0.62784615) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49470529) q[1];
sx q[1];
rz(-2.8637619) q[1];
sx q[1];
rz(-1.7956514) q[1];
x q[2];
rz(-2.4423994) q[3];
sx q[3];
rz(-2.1302038) q[3];
sx q[3];
rz(-1.4143816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.10852854) q[2];
sx q[2];
rz(-1.4427002) q[2];
sx q[2];
rz(2.4326883) q[2];
rz(2.3809643) q[3];
sx q[3];
rz(-1.1516738) q[3];
sx q[3];
rz(2.2061548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(3.0660504) q[0];
sx q[0];
rz(-1.1533371) q[0];
sx q[0];
rz(-0.72108889) q[0];
rz(-2.1636294) q[1];
sx q[1];
rz(-2.5826192) q[1];
sx q[1];
rz(1.5616547) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6128591) q[0];
sx q[0];
rz(-1.6472784) q[0];
sx q[0];
rz(-1.6441109) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97550895) q[2];
sx q[2];
rz(-1.3689318) q[2];
sx q[2];
rz(1.6109811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4750907) q[1];
sx q[1];
rz(-1.9618841) q[1];
sx q[1];
rz(0.60152454) q[1];
rz(-0.9004325) q[3];
sx q[3];
rz(-1.8547025) q[3];
sx q[3];
rz(2.8395124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28479031) q[2];
sx q[2];
rz(-0.52758104) q[2];
sx q[2];
rz(1.6214726) q[2];
rz(2.704845) q[3];
sx q[3];
rz(-2.2468061) q[3];
sx q[3];
rz(2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0590416) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(2.3178597) q[0];
rz(2.0027022) q[1];
sx q[1];
rz(-1.7887807) q[1];
sx q[1];
rz(1.9653856) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7438854) q[0];
sx q[0];
rz(-1.526014) q[0];
sx q[0];
rz(1.0683879) q[0];
rz(-pi) q[1];
rz(0.55032879) q[2];
sx q[2];
rz(-2.4541353) q[2];
sx q[2];
rz(-2.2000809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.215428) q[1];
sx q[1];
rz(-1.3326982) q[1];
sx q[1];
rz(1.5808616) q[1];
rz(-pi) q[2];
rz(2.5504712) q[3];
sx q[3];
rz(-1.6038461) q[3];
sx q[3];
rz(-0.37632468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99522432) q[2];
sx q[2];
rz(-0.58688846) q[2];
sx q[2];
rz(0.5874908) q[2];
rz(2.7557709) q[3];
sx q[3];
rz(-2.2930175) q[3];
sx q[3];
rz(-0.90252701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226444) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(2.5318085) q[0];
rz(1.3439641) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(-2.9885898) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333819) q[0];
sx q[0];
rz(-2.675288) q[0];
sx q[0];
rz(-3.0979041) q[0];
rz(3.1315003) q[2];
sx q[2];
rz(-0.79525012) q[2];
sx q[2];
rz(-0.45999664) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3045522) q[1];
sx q[1];
rz(-2.4992832) q[1];
sx q[1];
rz(-1.539597) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47346799) q[3];
sx q[3];
rz(-2.2064379) q[3];
sx q[3];
rz(1.9994232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0387705) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(0.62225303) q[2];
rz(-1.9410939) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(0.57435575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.95933652) q[0];
sx q[0];
rz(-3.0756364) q[0];
sx q[0];
rz(2.7863853) q[0];
rz(2.0390873) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(-2.4818518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099946215) q[0];
sx q[0];
rz(-2.2157241) q[0];
sx q[0];
rz(1.3553331) q[0];
rz(-pi) q[1];
rz(0.12730403) q[2];
sx q[2];
rz(-1.9812756) q[2];
sx q[2];
rz(0.32038996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5774571) q[1];
sx q[1];
rz(-0.45960765) q[1];
sx q[1];
rz(2.3338302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5054613) q[3];
sx q[3];
rz(-2.4628277) q[3];
sx q[3];
rz(-2.9836054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3416662) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(-0.1989092) q[2];
rz(0.79814664) q[3];
sx q[3];
rz(-1.1856368) q[3];
sx q[3];
rz(1.859349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2220919) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(-2.3566698) q[1];
sx q[1];
rz(-0.85826086) q[1];
sx q[1];
rz(0.8716743) q[1];
rz(-1.9368108) q[2];
sx q[2];
rz(-0.13199619) q[2];
sx q[2];
rz(0.97313626) q[2];
rz(2.8002515) q[3];
sx q[3];
rz(-0.69820709) q[3];
sx q[3];
rz(0.38521614) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
