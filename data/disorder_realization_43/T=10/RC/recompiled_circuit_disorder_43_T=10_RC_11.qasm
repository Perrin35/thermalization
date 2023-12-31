OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(-2.1656353) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358409) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(2.7266154) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33306723) q[2];
sx q[2];
rz(-2.1085848) q[2];
sx q[2];
rz(1.2531467) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.098903499) q[1];
sx q[1];
rz(-2.1747723) q[1];
sx q[1];
rz(0.16278111) q[1];
x q[2];
rz(-0.39335143) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(-2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7261937) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(-2.3577918) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.6911401) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(0.006342412) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3865249) q[0];
sx q[0];
rz(-1.7641983) q[0];
sx q[0];
rz(-0.13271876) q[0];
rz(0.96887529) q[2];
sx q[2];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.7768163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7087242) q[1];
sx q[1];
rz(-1.2607062) q[1];
sx q[1];
rz(-1.6354145) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3199602) q[3];
sx q[3];
rz(-2.1187966) q[3];
sx q[3];
rz(-2.8743924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94770849) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(0.095741622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.028713) q[0];
sx q[0];
rz(-2.9071147) q[0];
sx q[0];
rz(0.85588355) q[0];
rz(-2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(-0.72088748) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9853471) q[1];
sx q[1];
rz(-2.3932082) q[1];
sx q[1];
rz(-2.2845539) q[1];
x q[2];
rz(-2.5903969) q[3];
sx q[3];
rz(-1.6008198) q[3];
sx q[3];
rz(-0.72202819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.9018383) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-1.0338763) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(-2.1529249) q[0];
rz(-pi) q[1];
rz(-2.1448574) q[2];
sx q[2];
rz(-1.1166995) q[2];
sx q[2];
rz(-3.0388289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2305206) q[1];
sx q[1];
rz(-1.6358747) q[1];
sx q[1];
rz(1.3793524) q[1];
rz(-1.2437808) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(-0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(-0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-0.69798654) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(2.72686) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-0.57410747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78050437) q[0];
sx q[0];
rz(-1.6988806) q[0];
sx q[0];
rz(2.1696521) q[0];
x q[1];
rz(2.4262869) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(0.41358435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2327323) q[1];
sx q[1];
rz(-0.4545916) q[1];
sx q[1];
rz(0.67635398) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3247213) q[3];
sx q[3];
rz(-0.68021357) q[3];
sx q[3];
rz(-2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2888912) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067327499) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-0.69818991) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(1.408512) q[0];
rz(0.90768355) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(1.6931319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1440891) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-0.5703339) q[3];
sx q[3];
rz(-1.3864281) q[3];
sx q[3];
rz(-2.1031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7531062) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(-1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56918615) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(2.3715473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55243385) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(-0.60683672) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7355843) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(1.5232435) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8416482) q[1];
sx q[1];
rz(-0.83924676) q[1];
sx q[1];
rz(-2.4323835) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2392063) q[3];
sx q[3];
rz(-2.050403) q[3];
sx q[3];
rz(-0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-0.55505836) q[3];
sx q[3];
rz(-0.59593433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-0.24169895) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-0.48854488) q[1];
sx q[1];
rz(2.7752005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241656) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(-0.92745552) q[0];
x q[1];
rz(2.5744152) q[2];
sx q[2];
rz(-2.4293373) q[2];
sx q[2];
rz(0.30356193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.061325039) q[1];
sx q[1];
rz(-0.72239164) q[1];
sx q[1];
rz(-0.25218833) q[1];
x q[2];
rz(-1.0189692) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(0.036389694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-2.8477342) q[2];
rz(0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(-0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095187) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7991379) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(0.40823437) q[0];
rz(0.4653761) q[2];
sx q[2];
rz(-1.1306136) q[2];
sx q[2];
rz(-0.79681764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3190805) q[1];
sx q[1];
rz(-1.984239) q[1];
sx q[1];
rz(-0.49088571) q[1];
rz(-pi) q[2];
rz(-1.0320372) q[3];
sx q[3];
rz(-1.7146535) q[3];
sx q[3];
rz(2.0335576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.066594921) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(2.1203314) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-2.5773541) q[3];
sx q[3];
rz(-2.5436201) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026697712) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(0.23884493) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(-1.7609319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1191694) q[0];
sx q[0];
rz(-1.701231) q[0];
sx q[0];
rz(-0.93956691) q[0];
rz(-pi) q[1];
rz(2.0062749) q[2];
sx q[2];
rz(-1.5473817) q[2];
sx q[2];
rz(1.5362816) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7271991) q[1];
sx q[1];
rz(-2.6135923) q[1];
sx q[1];
rz(0.90331932) q[1];
x q[2];
rz(2.9525083) q[3];
sx q[3];
rz(-1.7366647) q[3];
sx q[3];
rz(-1.7961111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(-2.805368) q[2];
rz(-1.1296889) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(2.1375785) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-2.3168646) q[2];
sx q[2];
rz(-1.3168954) q[2];
sx q[2];
rz(-1.6864824) q[2];
rz(-1.929677) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
