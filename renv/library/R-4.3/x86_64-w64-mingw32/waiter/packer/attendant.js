!function(e,t){"object"==typeof exports&&"object"==typeof module?module.exports=t(require("Shiny"),require("jQuery")):"function"==typeof define&&define.amd?define(["Shiny","jQuery"],t):"object"==typeof exports?exports.attendant=t(require("Shiny"),require("jQuery")):e.attendant=t(e.Shiny,e.jQuery)}(self,(function(e,t){return(()=>{"use strict";var r={230:t=>{t.exports=e},145:e=>{e.exports=t}},a={};function n(e){var t=a[e];if(void 0!==t)return t.exports;var o=a[e]={exports:{}};return r[e](o,o.exports,n),o.exports}n.n=e=>{var t=e&&e.__esModule?()=>e.default:()=>e;return n.d(t,{a:t}),t},n.d=(e,t)=>{for(var r in t)n.o(t,r)&&!n.o(e,r)&&Object.defineProperty(e,r,{enumerable:!0,get:t[r]})},n.o=(e,t)=>Object.prototype.hasOwnProperty.call(e,t),n.r=e=>{"undefined"!=typeof Symbol&&Symbol.toStringTag&&Object.defineProperty(e,Symbol.toStringTag,{value:"Module"}),Object.defineProperty(e,"__esModule",{value:!0})};var o={};return(()=>{n.r(o),n(230),n(145);let e=new Map;const t=e=>{let t=$(`#${e.id} .progress-bar`),n=a(e.id,e.value);$(`#${e.id}`).show(),t.attr("aria-valuenow",e.value).css("width",n+"%"),e.text&&t.html(e.text),e.hideOnEnd&&r(e.id)<=e.value&&$(`#${e.id}`).hide()},r=e=>{let t=$(`#${e} .progress-bar`).attr("aria-valuemax");return parseFloat(t)},a=(e,t)=>t/r(e)*100;Shiny.addCustomMessageHandler("attendant-set-min",(e=>{$(`#${e.id} .progress-bar`).attr("aria-valuemin",e.min)})),Shiny.addCustomMessageHandler("attendant-set-max",(e=>{$(`#${e.id} .progress-bar`).attr("aria-valuemax",e.max)})),Shiny.addCustomMessageHandler("attendant-set",t),Shiny.addCustomMessageHandler("attendant-done",(a=>{let n=e.get(a.id);null!=n&&clearInterval(n),a.value=r(a.id),t(a)})),Shiny.addCustomMessageHandler("attendant-auto",(t=>{let n=$(`#${t.id} .progress-bar`),o=r(t.id),d=parseFloat(n.attr("aria-valuenow")),i=setInterval((()=>{if(d>o)return;d+=1;let e=a(t.id,d);$(`#${t.id} .progress-bar`).attr("aria-valuenow",d).css("width",e+"%")}),t.ms);e.set(t.id,i)}))})(),o})()}));