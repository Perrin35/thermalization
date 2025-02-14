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
rz(1.8472449) q[0];
sx q[0];
rz(-1.6142774) q[0];
sx q[0];
rz(-0.18728988) q[0];
rz(-1.3297431) q[1];
sx q[1];
rz(-2.5951374) q[1];
sx q[1];
rz(1.3475013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2696323) q[0];
sx q[0];
rz(-2.573421) q[0];
sx q[0];
rz(-1.4251891) q[0];
x q[1];
rz(2.7653469) q[2];
sx q[2];
rz(-2.0587789) q[2];
sx q[2];
rz(0.03586344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0986024) q[1];
sx q[1];
rz(-2.2523551) q[1];
sx q[1];
rz(-3.1206162) q[1];
rz(-pi) q[2];
rz(-0.25661664) q[3];
sx q[3];
rz(-1.768075) q[3];
sx q[3];
rz(2.9159604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61907855) q[2];
sx q[2];
rz(-1.9587025) q[2];
sx q[2];
rz(2.8600233) q[2];
rz(1.7742026) q[3];
sx q[3];
rz(-1.2308246) q[3];
sx q[3];
rz(0.27845964) q[3];
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
rz(-0.53681579) q[0];
sx q[0];
rz(-0.29412687) q[0];
sx q[0];
rz(-1.7915223) q[0];
rz(3.0916832) q[1];
sx q[1];
rz(-0.3388181) q[1];
sx q[1];
rz(-1.2443589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8690714) q[0];
sx q[0];
rz(-1.6289113) q[0];
sx q[0];
rz(-2.7921667) q[0];
rz(0.18644615) q[2];
sx q[2];
rz(-2.1121587) q[2];
sx q[2];
rz(0.10554927) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0383953) q[1];
sx q[1];
rz(-1.8382676) q[1];
sx q[1];
rz(2.9385376) q[1];
x q[2];
rz(-0.54173754) q[3];
sx q[3];
rz(-1.1835464) q[3];
sx q[3];
rz(0.70808402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0656978) q[2];
sx q[2];
rz(-1.1827129) q[2];
sx q[2];
rz(-0.069570216) q[2];
rz(-0.73404297) q[3];
sx q[3];
rz(-0.78248048) q[3];
sx q[3];
rz(0.61843425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065539) q[0];
sx q[0];
rz(-2.9422308) q[0];
sx q[0];
rz(2.8850436) q[0];
rz(1.8007295) q[1];
sx q[1];
rz(-0.95756617) q[1];
sx q[1];
rz(-2.0975013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0437735) q[0];
sx q[0];
rz(-1.7390842) q[0];
sx q[0];
rz(-1.2531444) q[0];
rz(-pi) q[1];
rz(-0.759941) q[2];
sx q[2];
rz(-2.8360998) q[2];
sx q[2];
rz(1.0585275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0699745) q[1];
sx q[1];
rz(-2.721792) q[1];
sx q[1];
rz(-2.384925) q[1];
rz(0.14117853) q[3];
sx q[3];
rz(-1.0021082) q[3];
sx q[3];
rz(1.6077667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9112245) q[2];
sx q[2];
rz(-2.0679097) q[2];
sx q[2];
rz(2.4252841) q[2];
rz(-0.45799842) q[3];
sx q[3];
rz(-1.674998) q[3];
sx q[3];
rz(-0.23076335) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79391795) q[0];
sx q[0];
rz(-2.1383998) q[0];
sx q[0];
rz(0.66116655) q[0];
rz(-1.8874774) q[1];
sx q[1];
rz(-1.5690119) q[1];
sx q[1];
rz(-1.5428146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5970478) q[0];
sx q[0];
rz(-2.8327541) q[0];
sx q[0];
rz(0.95914532) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2229052) q[2];
sx q[2];
rz(-2.4714758) q[2];
sx q[2];
rz(-2.5958259) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49816541) q[1];
sx q[1];
rz(-1.6339059) q[1];
sx q[1];
rz(-0.7038415) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38971735) q[3];
sx q[3];
rz(-1.0709312) q[3];
sx q[3];
rz(-0.12382774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25518146) q[2];
sx q[2];
rz(-0.35021439) q[2];
sx q[2];
rz(-2.0070455) q[2];
rz(-0.55401951) q[3];
sx q[3];
rz(-1.3628682) q[3];
sx q[3];
rz(-2.5567283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2968813) q[0];
sx q[0];
rz(-0.0021332707) q[0];
sx q[0];
rz(-1.1312436) q[0];
rz(3.0834815) q[1];
sx q[1];
rz(-1.2429712) q[1];
sx q[1];
rz(-1.6910472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3076271) q[0];
sx q[0];
rz(-1.183032) q[0];
sx q[0];
rz(0.8967171) q[0];
x q[1];
rz(1.8671473) q[2];
sx q[2];
rz(-0.83186281) q[2];
sx q[2];
rz(2.8339164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1102396) q[1];
sx q[1];
rz(-1.3340557) q[1];
sx q[1];
rz(-2.8048022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7282894) q[3];
sx q[3];
rz(-2.163475) q[3];
sx q[3];
rz(2.1324219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45945534) q[2];
sx q[2];
rz(-2.5294999) q[2];
sx q[2];
rz(-1.7876145) q[2];
rz(2.5891384) q[3];
sx q[3];
rz(-2.2508636) q[3];
sx q[3];
rz(-2.6433105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246493) q[0];
sx q[0];
rz(-2.1857078) q[0];
sx q[0];
rz(1.0981052) q[0];
rz(1.8662628) q[1];
sx q[1];
rz(-1.2397436) q[1];
sx q[1];
rz(1.0414418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5811036) q[0];
sx q[0];
rz(-2.1745918) q[0];
sx q[0];
rz(0.93224705) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0225811) q[2];
sx q[2];
rz(-1.6023984) q[2];
sx q[2];
rz(-2.1181483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66646092) q[1];
sx q[1];
rz(-1.0589927) q[1];
sx q[1];
rz(0.11630897) q[1];
rz(-2.8547229) q[3];
sx q[3];
rz(-0.86599819) q[3];
sx q[3];
rz(-1.6861841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2991403) q[2];
sx q[2];
rz(-1.9879397) q[2];
sx q[2];
rz(-0.057272043) q[2];
rz(2.9412269) q[3];
sx q[3];
rz(-0.55557957) q[3];
sx q[3];
rz(2.949775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6054691) q[0];
sx q[0];
rz(-2.496688) q[0];
sx q[0];
rz(0.38240018) q[0];
rz(2.1740225) q[1];
sx q[1];
rz(-0.41125527) q[1];
sx q[1];
rz(1.254902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9318274) q[0];
sx q[0];
rz(-1.5146274) q[0];
sx q[0];
rz(0.61135354) q[0];
x q[1];
rz(2.228565) q[2];
sx q[2];
rz(-2.0951197) q[2];
sx q[2];
rz(-1.9425478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6204024) q[1];
sx q[1];
rz(-0.56439059) q[1];
sx q[1];
rz(-2.317313) q[1];
rz(-pi) q[2];
rz(2.8627235) q[3];
sx q[3];
rz(-2.3430716) q[3];
sx q[3];
rz(2.4912733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.81662285) q[2];
sx q[2];
rz(-0.35494706) q[2];
sx q[2];
rz(-0.85186446) q[2];
rz(-0.24168092) q[3];
sx q[3];
rz(-1.9342187) q[3];
sx q[3];
rz(2.2309301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8082751) q[0];
sx q[0];
rz(-2.9589544) q[0];
sx q[0];
rz(2.0727378) q[0];
rz(1.3775657) q[1];
sx q[1];
rz(-0.27605468) q[1];
sx q[1];
rz(2.4527803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65387621) q[0];
sx q[0];
rz(-2.6453995) q[0];
sx q[0];
rz(-0.74183334) q[0];
rz(0.60283349) q[2];
sx q[2];
rz(-1.4571462) q[2];
sx q[2];
rz(1.08504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88676597) q[1];
sx q[1];
rz(-2.6313562) q[1];
sx q[1];
rz(-0.6986104) q[1];
x q[2];
rz(2.6121871) q[3];
sx q[3];
rz(-1.7200618) q[3];
sx q[3];
rz(2.4594027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0869202) q[2];
sx q[2];
rz(-2.4189147) q[2];
sx q[2];
rz(-1.5642536) q[2];
rz(2.4857322) q[3];
sx q[3];
rz(-1.4411996) q[3];
sx q[3];
rz(-0.76914579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5250788) q[0];
sx q[0];
rz(-2.109313) q[0];
sx q[0];
rz(0.14933625) q[0];
rz(-2.0556045) q[1];
sx q[1];
rz(-0.95567742) q[1];
sx q[1];
rz(0.48233262) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93123589) q[0];
sx q[0];
rz(-1.4853059) q[0];
sx q[0];
rz(-2.0443127) q[0];
rz(3.1343757) q[2];
sx q[2];
rz(-1.8033779) q[2];
sx q[2];
rz(-0.54534334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9378523) q[1];
sx q[1];
rz(-2.3934869) q[1];
sx q[1];
rz(2.7341873) q[1];
x q[2];
rz(-2.7217322) q[3];
sx q[3];
rz(-1.4814988) q[3];
sx q[3];
rz(1.4858703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0444191) q[2];
sx q[2];
rz(-2.3046875) q[2];
sx q[2];
rz(-2.6824299) q[2];
rz(3.0432213) q[3];
sx q[3];
rz(-1.8561074) q[3];
sx q[3];
rz(1.2488825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2761053) q[0];
sx q[0];
rz(-1.2609755) q[0];
sx q[0];
rz(3.0011445) q[0];
rz(1.9239931) q[1];
sx q[1];
rz(-2.0001037) q[1];
sx q[1];
rz(-1.5906895) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80181304) q[0];
sx q[0];
rz(-2.411541) q[0];
sx q[0];
rz(3.0193437) q[0];
rz(-pi) q[1];
rz(1.7459007) q[2];
sx q[2];
rz(-1.7410472) q[2];
sx q[2];
rz(-2.424602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1223476) q[1];
sx q[1];
rz(-0.68210318) q[1];
sx q[1];
rz(-2.7191799) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.012926558) q[3];
sx q[3];
rz(-2.2416246) q[3];
sx q[3];
rz(2.304043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6259674) q[2];
sx q[2];
rz(-2.1145861) q[2];
sx q[2];
rz(-1.94858) q[2];
rz(0.48202816) q[3];
sx q[3];
rz(-2.7218282) q[3];
sx q[3];
rz(-2.1217864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1358418) q[0];
sx q[0];
rz(-2.0477722) q[0];
sx q[0];
rz(2.3070106) q[0];
rz(0.76345481) q[1];
sx q[1];
rz(-1.3044985) q[1];
sx q[1];
rz(0.34674092) q[1];
rz(-1.687526) q[2];
sx q[2];
rz(-1.5925831) q[2];
sx q[2];
rz(1.5769373) q[2];
rz(-0.59178036) q[3];
sx q[3];
rz(-0.68282489) q[3];
sx q[3];
rz(0.015431701) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
