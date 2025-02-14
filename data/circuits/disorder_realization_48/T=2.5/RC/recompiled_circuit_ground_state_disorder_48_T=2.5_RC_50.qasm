OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(-0.40325525) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75664038) q[0];
sx q[0];
rz(-2.1813574) q[0];
sx q[0];
rz(-1.8288307) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4202021) q[2];
sx q[2];
rz(-0.70993844) q[2];
sx q[2];
rz(0.89872724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.510281) q[1];
sx q[1];
rz(-1.3461539) q[1];
sx q[1];
rz(2.2366877) q[1];
x q[2];
rz(0.81108629) q[3];
sx q[3];
rz(-2.1032277) q[3];
sx q[3];
rz(-2.6528751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.69748059) q[2];
sx q[2];
rz(-2.8287502) q[2];
sx q[2];
rz(-0.035813896) q[2];
rz(-0.73016417) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(-0.20619503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25074729) q[0];
sx q[0];
rz(-2.146281) q[0];
sx q[0];
rz(-0.19749755) q[0];
rz(0.47710553) q[1];
sx q[1];
rz(-2.4767866) q[1];
sx q[1];
rz(-0.8055996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2458513) q[0];
sx q[0];
rz(-0.68994683) q[0];
sx q[0];
rz(0.2941546) q[0];
rz(-pi) q[1];
rz(0.38061541) q[2];
sx q[2];
rz(-1.632649) q[2];
sx q[2];
rz(0.46910367) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.14368421) q[1];
sx q[1];
rz(-0.98947462) q[1];
sx q[1];
rz(0.79400009) q[1];
x q[2];
rz(-0.15112215) q[3];
sx q[3];
rz(-1.7100705) q[3];
sx q[3];
rz(2.3230011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3104559) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(1.7384701) q[2];
rz(-0.63203114) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(-0.12701756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3056575) q[0];
sx q[0];
rz(-0.62863612) q[0];
sx q[0];
rz(2.054731) q[0];
rz(-0.19733363) q[1];
sx q[1];
rz(-1.6355762) q[1];
sx q[1];
rz(-2.8089583) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98868552) q[0];
sx q[0];
rz(-1.1411785) q[0];
sx q[0];
rz(0.15362413) q[0];
rz(-pi) q[1];
rz(3.0789571) q[2];
sx q[2];
rz(-2.3047559) q[2];
sx q[2];
rz(0.96657414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3217993) q[1];
sx q[1];
rz(-1.395257) q[1];
sx q[1];
rz(1.3513397) q[1];
rz(-0.72686355) q[3];
sx q[3];
rz(-1.3976025) q[3];
sx q[3];
rz(-1.3722591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7670333) q[2];
sx q[2];
rz(-1.5908396) q[2];
sx q[2];
rz(-1.4683051) q[2];
rz(0.0091920216) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(0.79206842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.62995768) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(1.6988423) q[0];
rz(2.1171782) q[1];
sx q[1];
rz(-1.2316278) q[1];
sx q[1];
rz(-2.3146497) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073155135) q[0];
sx q[0];
rz(-0.19864635) q[0];
sx q[0];
rz(1.9563) q[0];
rz(-pi) q[1];
rz(-0.88281544) q[2];
sx q[2];
rz(-1.6903631) q[2];
sx q[2];
rz(-1.9012698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1093028) q[1];
sx q[1];
rz(-0.48498165) q[1];
sx q[1];
rz(0.46320199) q[1];
rz(-pi) q[2];
rz(-0.12110658) q[3];
sx q[3];
rz(-0.89205974) q[3];
sx q[3];
rz(-0.27287441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3247165) q[2];
sx q[2];
rz(-0.78410316) q[2];
sx q[2];
rz(1.2775705) q[2];
rz(-2.4534524) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(-2.1791606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0012896) q[0];
sx q[0];
rz(-1.4411417) q[0];
sx q[0];
rz(0.21743123) q[0];
rz(0.28542074) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(2.6434456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1284243) q[0];
sx q[0];
rz(-1.2449055) q[0];
sx q[0];
rz(-0.89114916) q[0];
rz(2.010538) q[2];
sx q[2];
rz(-0.56098191) q[2];
sx q[2];
rz(1.2351241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0302311) q[1];
sx q[1];
rz(-1.3627137) q[1];
sx q[1];
rz(-1.0743121) q[1];
rz(-0.50156939) q[3];
sx q[3];
rz(-2.6488228) q[3];
sx q[3];
rz(-1.2214965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.51776) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(-2.9916054) q[2];
rz(0.013966694) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(2.6278031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(-2.3288222) q[0];
rz(-1.7236408) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(2.0779804) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3115754) q[0];
sx q[0];
rz(-0.63236134) q[0];
sx q[0];
rz(2.564179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16367775) q[2];
sx q[2];
rz(-1.5541758) q[2];
sx q[2];
rz(-1.7134242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2460766) q[1];
sx q[1];
rz(-0.59883307) q[1];
sx q[1];
rz(3.1213698) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.205414) q[3];
sx q[3];
rz(-2.8293316) q[3];
sx q[3];
rz(-2.0978417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.873988) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(-1.7248636) q[2];
rz(-3.0806165) q[3];
sx q[3];
rz(-2.2305326) q[3];
sx q[3];
rz(-1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-3.022497) q[0];
rz(-2.6630317) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(0.68971577) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0478435) q[0];
sx q[0];
rz(-2.1458573) q[0];
sx q[0];
rz(-2.439365) q[0];
x q[1];
rz(1.941626) q[2];
sx q[2];
rz(-1.1111819) q[2];
sx q[2];
rz(2.4209765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3773681) q[1];
sx q[1];
rz(-2.2768094) q[1];
sx q[1];
rz(3.0824667) q[1];
rz(-0.85817091) q[3];
sx q[3];
rz(-1.458805) q[3];
sx q[3];
rz(1.5369161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0488284) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(-1.0678585) q[2];
rz(-2.974406) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(0.13944496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(-0.90840489) q[0];
rz(2.1482229) q[1];
sx q[1];
rz(-2.9790331) q[1];
sx q[1];
rz(2.2672674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4607166) q[0];
sx q[0];
rz(-1.5534288) q[0];
sx q[0];
rz(-0.13219035) q[0];
rz(-pi) q[1];
rz(-1.4972009) q[2];
sx q[2];
rz(-0.25840595) q[2];
sx q[2];
rz(-1.1272573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.43046865) q[1];
sx q[1];
rz(-1.6301155) q[1];
sx q[1];
rz(1.435168) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0197952) q[3];
sx q[3];
rz(-1.8809109) q[3];
sx q[3];
rz(1.9581025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5270762) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(-0.85421526) q[2];
rz(-2.6863875) q[3];
sx q[3];
rz(-0.56536094) q[3];
sx q[3];
rz(0.2555041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.47827569) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(1.664337) q[0];
rz(1.5669589) q[1];
sx q[1];
rz(-0.68395558) q[1];
sx q[1];
rz(-2.4050567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111209) q[0];
sx q[0];
rz(-2.6101042) q[0];
sx q[0];
rz(2.0474252) q[0];
rz(0.38072724) q[2];
sx q[2];
rz(-1.4743311) q[2];
sx q[2];
rz(-2.7166933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.29727259) q[1];
sx q[1];
rz(-1.2579009) q[1];
sx q[1];
rz(-1.827153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6227342) q[3];
sx q[3];
rz(-0.59762663) q[3];
sx q[3];
rz(-2.290987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3972724) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.6142023) q[2];
rz(-0.39198908) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.117332) q[0];
sx q[0];
rz(-1.6680822) q[0];
sx q[0];
rz(-0.20183739) q[0];
rz(2.9182538) q[1];
sx q[1];
rz(-2.6170862) q[1];
sx q[1];
rz(0.39252678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4215764) q[0];
sx q[0];
rz(-0.47552738) q[0];
sx q[0];
rz(-1.8029193) q[0];
x q[1];
rz(3.0891339) q[2];
sx q[2];
rz(-2.3067368) q[2];
sx q[2];
rz(-0.8142161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10822693) q[1];
sx q[1];
rz(-1.997233) q[1];
sx q[1];
rz(0.99797499) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48247561) q[3];
sx q[3];
rz(-1.2996718) q[3];
sx q[3];
rz(-1.1678054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7212123) q[2];
sx q[2];
rz(-0.98782295) q[2];
sx q[2];
rz(-0.58050275) q[2];
rz(2.765559) q[3];
sx q[3];
rz(-0.97085634) q[3];
sx q[3];
rz(-1.6002801) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15039438) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(-1.4906384) q[1];
sx q[1];
rz(-1.8883659) q[1];
sx q[1];
rz(-2.2813003) q[1];
rz(-0.80581325) q[2];
sx q[2];
rz(-1.0849107) q[2];
sx q[2];
rz(2.8203865) q[2];
rz(3.0720465) q[3];
sx q[3];
rz(-2.8290717) q[3];
sx q[3];
rz(1.1858218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
