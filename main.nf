#!/usr/bin/env nextflow

/*Pipeline de prueba para  nextflow */

/* Activar modo DSL2*/
nextflow.enable.dsl=2

params.str = 'Hello world!'

str_ch = Channel.from(params.str)

/*
 * Import modules
 */
include {
  splitLetters;
  convertToUpper } from './nf_modules/modules.nf'

/*
 * main pipeline logic
 */

workflow  {
  letters_ch = splitLetters(str_ch)
  upper_ch = convertToUpper(letters_ch.flatten())
  upper_ch.view{ it.trim()}
}
